#include<random>
#include "omp.h"
#include<algorithm>

// Parallelize using OMP

class grad_desc
{
    std::vector<int>                lbl;        // vector of labels (-1 or 1)
    std::vector< std::vector<int> > gnz_loc;    // matrix of locations of nonzero elements
    std::vector<int>                lnz_loc;    // vector of locations IN ONE ROW of nonzero elements
    std::vector< std::vector<double> > gnz_val; // matrix of values of nonzero elements
    std::vector<double>             lnz_val;    // vector of values IN ONE ROW of nonzero elements
    std::vector<double>             beta;       // concatenated solution vector
    std::vector<double>             beta_old;   // concatenated solution vector previous time step
    int                             lab;        // not sure?
    double                          mu;         // regularizer
    double                          N_d;        // number of data points (num rows)
    double                          N_f;        // number of data features (num cols)

public:
    grad_desc(std::string );
    double compute_global_f();
    double compute_global_f_me();
    double compute_beta_dot_x(int, int);
    double compute_beta_old_dot_x(int,int);
    double compute_sigma_exp_beta_dot_x(int);
    double compute_sigma_exp_beta_old_dot_x(int);
    std::vector<double> compute_grad_f_local_gd(int,bool);
    std::vector<double> compute_grad_f_local_sgd(int,bool);
    std::vector<double> compute_grad_f_global(bool);
    double compute_l2_norm_beta();
    int compute_support_beta();
    void print_beta();
    double choose_step_size(int);
    void update_beta_old();
    void perform_gd(int );
    void perform_one_step_gd(double);
    void perform_one_step_sgd(double);
    void perform_one_step_svrgd(double,std::vector<double>);
    void perform_parallel_sgd(int );
    void perform_parallel_svrgd(int);
    int random_data_point();
};

// initializes object and reads data from file
grad_desc::grad_desc(std::string data_fname)
{
    std::ifstream infile;
    infile.open(data_fname.c_str());
    std::string line,data,token;
    int count=0;
    mu=1e-8; 

    while(std::getline(infile,line))
    {
        std::istringstream iss(line);
	iss>>data;
        lab=std::stoi(data);
	lbl.push_back(lab);
	while(iss>>data)
	{
	    std::istringstream lss(data);
            while(std::getline(lss,token,':'))
	    {
	       if(count%2==0)
	       {
	           lnz_loc.push_back(std::stoi(token)-1);
	       }
	       else
	       {
	           lnz_val.push_back(std::stod(token));
	       }
	       count++;
	    }
	}
	count=0;
	gnz_loc.push_back(lnz_loc);
	gnz_val.push_back(lnz_val);
	lnz_loc.clear();
	lnz_val.clear();
    }

    /*for(int i=0;i<lbl.size();i++)
    {
        std::cout<<"Label "<<lbl[i]<<std::endl;
	for(int j=0;j<gnz_loc[i].size();j++)
	{
	    std::cout<<gnz_loc[i][j]<<":"<<gnz_val[i][j]<<" ";
	}
	std::cout<<std::endl;
    }*/
    N_f=47236;
    beta.resize(2*N_f);
    beta_old.resize(2*N_f);
    N_d=lbl.size();
}

// compute L2 norm of the full beta vector
double grad_desc::compute_l2_norm_beta()
{
    double res=0.0;
    for(int i=0;i<beta.size();i++)
    {
	res=res + beta[i]*beta[i]; 
    }
    return res;
}

// compute the support of the full beta vector
int grad_desc::compute_support_beta(){
    int res=0;
    for(int i=0;i<beta.size();i++)
    {
        if (beta[i]!=0.0) res++;
    }
    return res;
}

// print the full beta vector to terminal
void grad_desc::print_beta(){
    for(int i=0;i<beta.size();i++)
    {
        std::cout << "beta: " << beta[i] << std::endl;
    }
    return;
}

// compute the global value of f
// with regularizer included
double grad_desc::compute_global_f()
{
    double res=0.0;
    int idx;

    for(int i=0;i<N_d;i++)
    {
        if(lbl[i]==-1)
    	{
    	    idx=0;
    	}
    	else
    	{
    	    idx=1;
    	}
        res= res + (compute_beta_dot_x(idx,i) + log(compute_sigma_exp_beta_dot_x(i))); 
    }
    res=res/N_d;
    res=res+mu*compute_l2_norm_beta();
    return res;
}

// compute the global value of f
// WITHOUT the regularizer
double grad_desc::compute_global_f_me()
{
    double res=0.0;
    int idx;

    for(int i=0;i<N_d;i++)
    {
        if(lbl[i]==-1)
	{
	    idx=0;
	}
	else
	{
	    idx=1;
	}
        res= res + (compute_beta_dot_x(idx,i) + log(compute_sigma_exp_beta_dot_x(i))); 
    }
    //res=res + mu*compute_l2_norm_beta();
    res=res/N_d;
    return res;
}


double grad_desc::compute_beta_dot_x(int idx,int i)
{
    double res=0.0;
    for(int j=0;j<gnz_loc[i].size();j++)
    {
        res=res + beta[idx*N_f+gnz_loc[i][j]]*gnz_val[i][j];
    }
    return res;
}

void grad_desc::update_beta_old()
{
    // Copy the current beta and store it into beta_old

    std::copy(beta.begin(),beta.end(),beta_old.begin());
}

double grad_desc::compute_beta_old_dot_x(int idx,int i)
{
    double res=0.0;
    for(int j=0;j<gnz_loc[i].size();j++)
    {
        res=res + beta_old[idx*N_f+gnz_loc[i][j]]*gnz_val[i][j];
    }
    return res;
}

double grad_desc::compute_sigma_exp_beta_dot_x(int i)
{
    double res=0.0;
    res=res + exp(-1*compute_beta_dot_x(0,i)) + exp(-1*compute_beta_dot_x(1,i));
    return res;
}

double grad_desc::compute_sigma_exp_beta_old_dot_x(int i)
{
    double res=0.0;
    res=res + exp(-1*compute_beta_old_dot_x(0,i)) + exp(-1*compute_beta_old_dot_x(1,i));
    return res;
}

// compute gradient of local f (one data point)
// EXCLUDING the regularizer term
std::vector<double>  grad_desc::compute_grad_f_local_sgd(int i, bool BETA_FLAG)
{
    std::vector<double> grad;
    grad.resize(N_f*2);
    int base=0;
    double N=lbl.size();

    //Adding contribution to gradient from log term
    for(int j=0;j<gnz_loc[i].size();j++)
    {
    	if(BETA_FLAG == 1)
    	{
    	    grad[gnz_loc[i][j]]=-gnz_val[i][j]*exp(-compute_beta_dot_x(0,i))/(compute_sigma_exp_beta_dot_x(i));
    	    grad[gnz_loc[i][j]+N_f]=-gnz_val[i][j]*exp(-compute_beta_dot_x(1,i))/(compute_sigma_exp_beta_dot_x(i));
    	}
    	else
    	{

    	    grad[gnz_loc[i][j]]=-gnz_val[i][j]*exp(-compute_beta_old_dot_x(0,i))/(compute_sigma_exp_beta_old_dot_x(i));
    	    grad[gnz_loc[i][j]+N_f]=-gnz_val[i][j]*exp(-compute_beta_old_dot_x(1,i))/(compute_sigma_exp_beta_old_dot_x(i));
    	}
    }

    //Adding contribution to gradient from linear term
    if(lbl[i]==-1)
    {
         for(int j=0;j<gnz_loc[i].size();j++)
    	 {
    	     grad[gnz_loc[i][j]]=grad[gnz_loc[i][j]] + gnz_val[i][j];
    	 }
    }
    else
    {
        for(int j=0;j<gnz_loc[i].size();j++)
    	{
    	    grad[gnz_loc[i][j]+N_f]=grad[gnz_loc[i][j]+N_f] + gnz_val[i][j];
    	}
    }


    //Adding contribution from regularizer
   
    /*for(int j=0;j<grad.size();j++)
    {
        grad[j]=grad[j] + 2*(mu)*N_d*beta[j]; 
    }*/
    return grad;
    
}

// compute gradient of local f 
std::vector<double> grad_desc::compute_grad_f_local_gd(int i, bool BETA_FLAG)
{
    std::vector<double> grad;
    grad.resize(N_f*2);
    int base=0;
    double N=lbl.size();

    //Adding contribution to gradient from log term
    for(int j=0;j<gnz_loc[i].size();j++)
    {
    	if(BETA_FLAG == 1)
    	{
    	    grad[gnz_loc[i][j]]=-gnz_val[i][j]*exp(-compute_beta_dot_x(0,i))/(compute_sigma_exp_beta_dot_x(i));
    	    grad[gnz_loc[i][j]+N_f]=-gnz_val[i][j]*exp(-compute_beta_dot_x(1,i))/(compute_sigma_exp_beta_dot_x(i));
    	}
    	else
    	{
    	    grad[gnz_loc[i][j]]=-gnz_val[i][j]*exp(-compute_beta_old_dot_x(0,i))/(compute_sigma_exp_beta_old_dot_x(i));
    	    grad[gnz_loc[i][j]+N_f]=-gnz_val[i][j]*exp(-compute_beta_old_dot_x(1,i))/(compute_sigma_exp_beta_old_dot_x(i));
    	}
    }

    //Adding contribution to gradient from linear term
    if(lbl[i]==-1)
    {
         for(int j=0;j<gnz_loc[i].size();j++)
	 {
	     grad[gnz_loc[i][j]]=grad[gnz_loc[i][j]] + gnz_val[i][j];
	 }
    }
    else
    {
        for(int j=0;j<gnz_loc[i].size();j++)
	{
	    grad[gnz_loc[i][j]+N_f]=grad[gnz_loc[i][j]+N_f] + gnz_val[i][j];
	}
    }
    
    //Adding contribution from regularizer
   
    return grad;
    
}

// compute gradient of global f function
// equivalent to the sum of all gradients
std::vector<double> grad_desc::compute_grad_f_global(bool BETA_FLAG)
{
   std::vector<double> grad,grad_l;
   grad.resize(N_f*2);
   double N=lbl.size();

   #pragma omp parallel for default(shared) private(grad_l)
   for(int i=0;i<lbl.size();i++)
   {
        // if (omp_get_thread_num() == 0){
        //     std::cout << "i: " << i << "/" << lbl.size() << std::endl;
        // }
       grad_l=compute_grad_f_local_gd(i,BETA_FLAG);
       for(int j=0;j<grad.size();j++)
       {
           grad[j]=grad[j]+grad_l[j]/N_d;
       }
   }

   #pragma omp parallel for default(shared)
   for(int i=0;i<beta.size();i++)
   {
        grad[i]=grad[i] + 2*mu*beta[i];
   }

   /*for(int j=0;j<grad.size();j++)
   {
       if(grad[j]!=0)
       {
           std::cout<<grad[j]<<std::endl;
       }
   }*/
    
   return grad; 
}

/*
std::vector<double> grad_desc::compute_grad_f_global(bool USE_BETA_OLD)
{
   std::vector<double> grad,grad_l;
   grad.resize(N_f*2);
   double N=lbl.size();
   for(int i=0;i<lbl.size();i++)
   {
       grad_l=compute_grad_f_local_gd(i,USE_BETA_OLD);
       for(int j=0;j<grad.size();j++)
       {
           grad[j]=grad[j]+grad_l[j]/N_d;
       }
   }
   for(int i=0;i<beta.size();i++)
   {
        grad[i]=grad[i] + 2*mu*beta[i];
   }

   for(int j=0;j<grad.size();j++)
   {
       if(grad[j]!=0)
       {
           std::cout<<grad[j]<<std::endl;
       }
   }
    
   return grad; 
}
*/

void grad_desc::perform_one_step_gd(double step_size)
{
    std::vector<double> grad;
    std::cout<<"--------------One step GD-------------------"<<std::endl;
    std::cout<<"gd_error "<<compute_global_f()<<std::endl;

    bool USE_BETA_CURRENT = 1;

    grad=compute_grad_f_global(USE_BETA_CURRENT);

    for(int i=0;i<beta.size();i++)
    {
        beta[i]=beta[i]  - step_size*grad[i];
    }
    //std::cout<<"Value after one step of gd "<<compute_global_f()<<std::endl;
    std::cout<<"--------------One step GD done---------------"<<std::endl;
    std::cout<<std::endl;
}

int grad_desc::random_data_point()
{
    // Choosing random index for calculating the gradient term for sgd

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(0,N_d-1);
    return dis(gen);
}
void grad_desc::perform_one_step_sgd(double step_size)
{
    //----- Move random number generation to a function ------//

    std::vector<double> grad;
    int rand_pt = random_data_point();

    //std::random_device rd;
    //std::mt19937 gen(rd());
    //std::uniform_int_distribution<> dis(0,N_d-1);
    //std::vector<double> grad;
    //std::cout<<"--------------One step SGD-----------------"<<std::endl;
    //std::cout<<"Random data point "<<dis(gen)<<std::endl;

    bool USE_BETA_CURRENT = 1;		

    grad=compute_grad_f_local_sgd(rand_pt,USE_BETA_CURRENT);

    for(int i=0;i<beta.size();i++)
    {
        grad[i]=grad[i] + 2*mu*beta[i];
    }
    for(int i=0;i<beta.size();i++)
    {
        beta[i]=beta[i]  - step_size*grad[i];
    }

    //std::cout<<"--------------One step SGD performed-------"<<std::endl;
}

void grad_desc::perform_gd(int steps)
{
    int gd_steps=0;
    double step_size=100.0;
    while(gd_steps<steps)
    {
        std::cout<<"gd_step"<<gd_steps<<std::endl;
        perform_one_step_gd(step_size);
	gd_steps++;
    }
}


double grad_desc::choose_step_size(int gd_steps)
{

	double step_size;

	step_size = 5;

	if(gd_steps>500){
	step_size=2.5;
	}	
	if(gd_steps>5000){
	step_size=1.0;
	}
	if(gd_steps>200000){
	step_size=0.5;
	}
	if(gd_steps>300000){
	step_size=0.05;
	}
	if(gd_steps>400000){
	step_size=0.01;
	}

	return step_size;

}


void grad_desc::perform_parallel_sgd(int steps)
{
    int gd_steps=0;
    double step_size=5;
    double temp,opt_err;
    std::cout<<"Value before SGD "<<compute_global_f()<<std::endl;

    int num_threads;
    double t_begin,t_end,elapsed,tid;

    t_begin = omp_get_wtime();

    std::cout << "grepforme:, Iteration, Elapsed_Time, SGD_Error, Opt_Error, Norm_Beta, Sup_Beta" << std::endl;

    #pragma omp parallel private(tid)				// Spawns a group of threads
    {		 
    	tid = omp_get_thread_num();
	num_threads = omp_get_num_threads();

    # pragma omp for schedule(static)    					// divided the work amongst the spawned threads    	    
	for(gd_steps=0;gd_steps<steps;gd_steps++)
	{
	    //std::cout<<"sgd_step "<<gd_steps<<std::endl;
	    
	    perform_one_step_sgd(step_size);
	    step_size = choose_step_size(gd_steps);
	   

        //temp=compute_global_f_me();
	//opt_err=temp + mu*(compute_l2_norm_beta())/N_d;
        //std::cout<<compute_l2_norm_beta()<<std::endl;

	    if(tid==0){
//		std::cout << "gd_steps : \t" << gd_steps << std::endl;
	    if((gd_steps%1)==0){
		    double t_current,t_elapsed;

		    t_current = omp_get_wtime();
		    t_elapsed = t_current - t_begin;

		    temp=compute_global_f_me();
		    opt_err=temp + mu*(compute_l2_norm_beta());
		    // std::cout<<compute_l2_norm_beta()<<std::endl;
		    // std::cout<<"Step number \t"<<gd_steps<<"\t sgd_error \t"<<temp<<"\t opt_error \t"<<opt_err<<std::endl;
		    // std::cout<<"Time Elapsed :\t"<<t_elapsed << "  seconds" << std::endl;
	       std::cout << "grepforme:, " << gd_steps << ", " << t_elapsed << ", " ;
           std::cout << temp << ", " << opt_err << ", " << compute_l2_norm_beta() << ", " << compute_support_beta() << std::endl;
        }
	    }
	} // End of Pragma omp for
    }//End of Pragma omp parallel

    t_end = omp_get_wtime();
    elapsed = t_end-t_begin;

    print_beta();

    std::cout<<"Value after Parallel SGD "<<compute_global_f()<<std::endl;
    std::cout<<"Total Time Taken for \t"<<steps<<"\t Parallel sgd_steps is \t " << elapsed << "  seconds" << std::endl; 
}

void grad_desc::perform_one_step_svrgd(double step_size, std::vector<double>grad_global)
{
    std::vector<double>grad_sgd;
    std::vector<double>grad_sgd_with_beta_old;
    std::vector<double>grad_svrgd;

    bool USE_BETA_OLD = 0;			// Flag which instructs the function to use the old beta
    bool USE_BETA_CURRENT = 1;			// Flag which instructs the function to use the current beta

    grad_sgd.resize(N_f*2);
    grad_sgd_with_beta_old.resize(N_f*2);
    grad_svrgd.resize(N_f*2);

    int rand_pt = random_data_point();
    //std::cout << "Random \t" << rand_pt << std::endl;
    grad_sgd = compute_grad_f_local_sgd(rand_pt,USE_BETA_CURRENT);		// Calculated at the current beta FLAG == 0
    grad_sgd_with_beta_old = compute_grad_f_local_sgd(rand_pt,USE_BETA_OLD);	// Calculated using the old beta FLAG == 1

    for(int i=0;i<beta.size();i++)
    {
    	grad_svrgd[i] = 1.0*(grad_sgd[i] - grad_sgd_with_beta_old[i] + grad_global[i]);
    	beta[i] = beta[i]  - step_size*grad_svrgd[i];
    }

}

void grad_desc::perform_parallel_svrgd(int steps)
{
    
    int global_gd_steps,sgd_steps,global_gd_interval;
    double step_size=5;
    double temp,opt_err;

    #pragma omp parallel
    {
    #pragma omp single
    {
    std::cout << "Using " << omp_get_num_threads() << " threads" << std::endl;    
    std::cout<<" Value before SVRGD "<<compute_global_f()<<std::endl;
    }

    }

    std::vector<double>grad_global;
    grad_global.resize(N_f*2);

    global_gd_interval = (4e-2)*steps;			// Interval between global full gradient computation
    std::cout << "Global GD Interval \t" << global_gd_interval << std::endl;

    bool USE_BETA_OLD = 0;				// Flag which instructs the function to use the old beta
    bool USE_BETA_CURRENT = 1;				// Flag which instructs the function to use the current beta

    // Outer global gradient loop

    int tid;
    double t_begin,t_end,elapsed;

    t_begin = omp_get_wtime();
    
    std::cout << "grepforme:, Iteration, Elapsed_Time, SGD_Error, Opt_Error, Norm_Beta, Sup_Beta" << std::endl;

	
    for(global_gd_steps=0;global_gd_steps<steps;global_gd_steps+=global_gd_interval)
    {

	#pragma omp parallel private(tid) shared(grad_global)
	{	
    
	tid = omp_get_thread_num();	
	
	// Only use one thread for the outerloop
	

                std::vector<double> grad_l;
                grad_global.resize(N_f*2);
                double N=lbl.size();

               #pragma omp for
               for(int i=0;i<lbl.size();i++)
               {
                    // if (omp_get_thread_num() == 0){
                    //     std::cout << "i: " << i << "/" << lbl.size() << std::endl;
                    // }
                   grad_l=compute_grad_f_local_gd(i,USE_BETA_CURRENT);
                   for(int j=0;j<grad_global.size();j++)
                   {
                       grad_global[j]=grad_global[j]+grad_l[j]/N_d;
                   }
               }

               #pragma omp for 
               for(int i=0;i<beta.size();i++)
               {
                    grad_global[i]=grad_global[i] + 2*mu*beta[i];
               }

// 	#pragma omp single
// 	{
// 	    grad_global = compute_grad_f_global(USE_BETA_CURRENT);
// 	    update_beta_old();
// //	    std::cout <<" Global Gradient Computed" << std::endl;
// 	}

	// Inner SGD loop using the outer global gradient
    
	#pragma omp for 
	    for(sgd_steps = 0;sgd_steps<global_gd_interval;sgd_steps++)
	    {
		step_size = choose_step_size(global_gd_steps+sgd_steps);
		perform_one_step_svrgd(step_size,grad_global);			// updates the beta (solution) using both global and local gradient

		if(tid==0){
		if(((global_gd_steps+sgd_steps)%1)==0)
		{
		    double t_current,t_elapsed;

		    t_current = omp_get_wtime();
		    t_elapsed = t_current - t_begin;

		    temp=compute_global_f_me();
		    opt_err=temp + mu*(compute_l2_norm_beta());
            std::cout << "grepforme:, " << global_gd_steps+sgd_steps << ", " << t_elapsed << ", " ;
            std::cout << temp << ", " << opt_err << ", " << compute_l2_norm_beta() << ", " << compute_support_beta() << std::endl;

		    // std::cout<<compute_l2_norm_beta()<<std::endl;
		    // std::cout<<"Step number \t"<<global_gd_steps+sgd_steps<<"\t svrgd_error \t"<<temp<<"\t opt_error \t"<<opt_err<<std::endl;
		    // std::cout<<"Time Elapsed :\t"<<t_elapsed << "  seconds" << std::endl;
		}
		}
	    }// End omp for loop

	}// End omp parallel loop
    }//end for loop

    t_end = omp_get_wtime();
    elapsed = t_end-t_begin;

    print_beta();

    std::cout<<"Value after Parallel SVRGD "<<compute_global_f()<<std::endl;
    std::cout<<"Total Time Taken for \t"<<steps<<"\t Parallel svrgd_steps is \t " << elapsed << "  seconds" << std::endl; 
    
}	   
	    
