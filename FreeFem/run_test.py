#!/usr/bin/python

from subprocess import call
import os,itertools

print("Running the Freefem simulation on various configurations")

# Create the parameters file 

def write_param_file(gstep,beta,Nmax,FE_type,plot_data) : 
	
	param_file = open("scripts/numerical_parameters.edp",'w')
	param_file.write("// Configuration file written by Python testing script\n");

	param_file.write("real gstep = " + str(gstep) + " ;\n");
	param_file.write("real beta = " + str(beta) + " ;\n");
	param_file.write("int Nmax = " + str(Nmax) + " ;\n");

	param_file.write("\n") ;

	param_file.write("fespace V2h(Th,"+FE_type+") ;\n");
	param_file.write("fespace V2in(Tin,"+FE_type+") ;\n");

	if(plot_data) :  
		param_file.write("bool plotData = true ;")  
	else : 
		param_file.write("bool plotData = false ;")

	param_file.close()

	# end of function

def run_freefem() :
	old_path = os.getcwd();
	os.chdir(old_path+"/scripts")

	call(["FreeFem++","body_2d_inverse_interp.edp"])

	os.chdir(old_path)

	# end of function

def plot_J() : 
	# Write gnuplot config file
	f = open("visualize.conf",'w')
	f.write("#GNUPLOT script written by Python testing script\n\n")
	f.write("set terminal postscript eps enhanced color\n")
	f.write("set autoscale\nset xlabel 'Iteration'\nset ylabel 'Cost functional'\nset style data lines\n")

	f.write("set output 'export/images/J'\n")
	f.write("plot 'export/J.txt' lw 3 lc 'blue' title 'Cost Functional'\n")

	f.write("set output 'export/images/tJ'\n")
	f.write("plot 'export/tJ.txt' lw 3 lc 'blue' title 'true Cost Functional'\n")

	f.close()

	# Run gnuplot and then remove the config file
	call(["gnuplot","visualize.conf"])
	call(["rm","visualize.conf"])

	# end of function

def p1_vs_p2(gstep,beta,Nmax) : 

	print "P1 vs P2 test "
	test_name = str(Nmax)+"_"+str(beta)+"_"+str(gstep)
	call(["mkdir","-p","export/trials/"+test_name])
	# Run the P1 case
	print "Running the P1 case... "
	write_param_file(gstep,beta,Nmax,"P1",False)
	run_freefem()
	call(["mv","export/J.txt","export/trials/"+test_name+"/J_P1.txt"])

	# Run the P2 case
	print "Running the P2 case... "
	write_param_file(gstep,beta,Nmax,"P2",False)
	run_freefem()
	call(["mv","export/J.txt","export/trials/"+test_name+"/J_P2.txt"])
	
	# Plot output data
	f = open("visualize.conf",'w')
	f.write("#GNUPLOT script written by Python testing script\n\n")
	f.write("set terminal postscript eps enhanced color\n")
	f.write("set autoscale\nset xlabel 'Iteration'\nset ylabel 'Cost functional'\nset style data lines\n")

	f.write("set output 'export/trials/"+test_name+"/P1vsP2'\n")
	f.write("plot 'export/trials/"+test_name+"/J_P1.txt' lw 3 lc 'blue' title 'Cost Functional P1' , ")
	f.write("'export/trials/"+test_name+"/J_P2.txt' lw 3 lc 'red' title 'Cost Functional P2 \n")

	f.close()

	# Run gnuplot and then remove the config file
	call(["gnuplot","visualize.conf"])
	call(["rm","visualize.conf"])

	# end of function



def run_multiple_trials(gsteps,betas,Nmaxs,FE_types,TEST_NAME) : 

	# Create output folder
	call(["mkdir","-p","export/trials/"+TEST_NAME])

	# Run Freefem script
	for (gstep,beta,Nmax,FE_type) in itertools.product(gsteps,betas,Nmaxs,FE_types) : 
		test_name = str(Nmax)+"_"+str(beta)+"_"+str(gstep)+"_"+FE_type
		# Write param file
		write_param_file(gstep,beta,Nmax,FE_type,False)
		run_freefem()
		call(["mv","export/J.txt","export/trials/"+TEST_NAME+"/"+test_name+"_J.txt"])
		call(["mv","export/tJ.txt","export/trials/"+TEST_NAME+"/"+test_name+"_tJ.txt"])
		
	# end of function

def plot_trials_result(gsteps,betas,Nmaxs,FE_types,TEST_NAME) : 
	# Create output folder
	call(["mkdir","-p","export/trials/"+TEST_NAME+"/images"])

	# Write visualize.conf file
	f = open("visualize.conf",'w')
	f.write("#GNUPLOT script written by Python testing script\n\n")
	f.write("set terminal postscript eps enhanced color\n")
	f.write("set autoscale\n")

	# Plot all Js
	f.write("set xlabel 'Iteration'\nset ylabel 'Cost functional'\nset style data lines\n")
	f.write("\n\n")
	for (gstep,beta,Nmax,FE_type) in itertools.product(gsteps,betas,Nmaxs,FE_types) : 
		test_name = str(Nmax)+"_"+str(beta)+"_"+str(gstep)+"_"+FE_type
		f.write("set output 'export/trials/"+TEST_NAME+"/images/"+test_name+"_J'\n")
		f.write("plot 'export/trials/"+TEST_NAME+"/"+test_name+"_J.txt' lw 3 lc 'blue' title 'Cost Functional'\n")




	f.close()

	# Run gnuplot and remove conf file
	call(["gnuplot","visualize.conf"])
	call(["rm","visualize.conf"])

	# end of function


if __name__ == '__main__' : 

	# p1_vs_p2(0.004,0.001,20)

	# run_multiple_trials([0.004, 0.008],[0.001],[5],["P1"],"Test_Python")

	plot_trials_result([0.004, 0.008],[0.001],[5],["P1"],"Test_Python")