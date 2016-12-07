#include<iostream>
#include<algorithm>
#include<fstream>
#include<vector>
#include<string>
#include<sstream>
#include<cmath>
#include "grad_class.h"

//double dot_prod(std::vector<std::vector<double> >beta,int N);


int main()
{

    
    grad_desc *gd=new grad_desc("rcv1_train.binary");
    std::cout<<"Finished processing data"<<std::endl;
    gd->perform_parallel_sgd(10000);
    //gd->perform_parallel_svrgd(10000);
    return 0;
}
