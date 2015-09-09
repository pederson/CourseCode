# short script to do a homework problems 5 and 6 in hw1

def euler(yn, h):
	return yn + h*yn

def scheme(yn, ynp1, h):
	return (2+3*h)*yn - ynp1
	#return yn + 2*h*ynp1

h = 0.025
t = range(int(1/h)+1)
for i in range(len(t)):
	t[i] = t[i]*h
y0 = 1

yn = y0
ynp1 = euler(yn, h)
print ("t: %1.3f    y: %1.3f" % (t[1], ynp1))
for i in range(len(t)-2):
	curi = i+2
	curval = scheme(yn, ynp1, h)

	# print output
	print ("t: %1.3f    y: %1.3f" % (t[curi], curval))
	# update values
	yn = ynp1
	ynp1 = curval
