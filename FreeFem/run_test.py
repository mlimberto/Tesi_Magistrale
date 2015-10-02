#!/usr/bin/python

from subprocess import call
import os

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


if __name__ == '__main__' : 

	write_param_file(0.004,0.001,15,"P1",False);

	run_freefem()