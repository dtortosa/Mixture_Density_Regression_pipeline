#!/usr/bin/env python3.8


#############################################################
########## SCRIPT FOR MODELLING IHS IN GEN WINDOWS ##########
#############################################################

#Model the association between iHS and genomic factors across genomic windows using a Mixture Density Regression.



###################
##### IMPORTS #####
###################

from __future__ import print_function, division
import pandas as pd
import numpy as np
import scipy as sp 
from scipy import optimize
from sklearn import preprocessing 
import matplotlib.pyplot as plt
import multiprocessing as mp
import subprocess
import glob
import os
from functools import reduce
import re
from scipy.special import expit



########################################
########## SET THE DIRECTORIES #########
########################################

#set the working directory
path_inside_container = '/'
path_starting_folder_raw = subprocess.run(['pwd'], stdout=subprocess.PIPE)
path_starting_folder_raw.stdout
encoding = 'utf-8'
path_starting_folder = str(path_starting_folder_raw.stdout, encoding)
path_starting_folder = path_starting_folder.rstrip()

#set the path for inputs
path_outside_container_inputs = [path_starting_folder, '/data']
path_outside_container_inputs = ''.join(path_outside_container_inputs)

#set the path for outputs
path_outside_container_outputs = [path_starting_folder, '/results']
path_outside_container_outputs = ''.join(path_outside_container_outputs)


#make some directories
os.system('mkdir -p ' + path_outside_container_outputs)
os.system('mkdir -p ' + path_outside_container_outputs + '/figures')
os.system('mkdir -p ' + path_outside_container_outputs + '/tables')



#############################################################
################## EXTRACT POPULATION NAMES #################
#############################################################

#First, join paths to create the complete path to the iHS data
path_to_ihs_data = ''.join([path_outside_container_inputs, '/mean_ihs_gene_windows/*'])

#list all the files in that folder (ihs value per population and chromosome) including the whole path
list_ihs_files_full_path = glob.glob(path_to_ihs_data) 

#remove the path to the file name
list_ihs_files_only_names = [os.path.basename(x) for x in list_ihs_files_full_path]

#extract the population name from all the files
list_pops = [x.split('_')[0] for x in list_ihs_files_only_names]

#select the unique cases
list_pops_unique = np.unique(list_pops)

#we have to change list_pops_unique, a nunmpy.darray, to a list, because if not the parallelization does not work.
list_pops_unique = list_pops_unique.tolist()


##make a dict for getting the correct pop abbreviations
pop_abbreviations = {"YRID": "YRI", "TSID": "TSI", "CEUD": "CEU", "CHBD": "CHB", "PELD": "PEL"}


##########################################################################
################## DEFINE THE FUNCTION TO BE PARALLELIZED ################
##########################################################################

#selected_pop = list_pops_unique[0]
def ihs_modelling_per_pop(selected_pop):

	############################
	##### PREPARE THE DATA #####
	############################
	
	##ihs values and the number of iHS data points
	#bind the full path to the iHS data for the selected population
	ihs_data_path = ''.join([path_outside_container_inputs, '/mean_ihs_gene_windows/', selected_pop, '_mean_ihs_gene_windows_final_v1.txt.gz'])
	ihs_n_points_path = ''.join([path_outside_container_inputs, '/mean_ihs_gene_windows/', selected_pop, '_n_ihs_gene_windows_final_v1.txt.gz'])
		#we also want the number of iHS data points per window as a predictor in the models.
	
	#load the iHS data for the selected population across all the window sizes
	ihs_data = pd.read_csv(ihs_data_path, compression='gzip', header=0, sep='\t', low_memory=False)
	ihs_n_points = pd.read_csv(ihs_n_points_path, compression='gzip', header=0, sep='\t', low_memory=False)
	
	
	##predictors
	#create a list with the final version of all confounding factors
	#First, join paths to create the complete path to the iHS data
	path_to_predictor_data = ''.join([path_outside_container_inputs, '/predictors_gene_windows/*.txt'])
	
	#list of the file in that folder (ihs value per population and chromosome) including the whole path
	list_predictors_files_full_path = glob.glob(path_to_predictor_data) 
	
	#save the name of the predictors separately. These are the names of the list of dataframes that we are going to create.
	list_predictors_files_only_names = [os.path.basename(x) for x in list_predictors_files_full_path]
	
	#load all the dataframes for which we have path within the list
	list_df_predictors = [pd.read_csv(x, header=0, sep='\t', low_memory=False) for x in list_predictors_files_full_path] 
	

	
	##############################################################
	################## OPEN LOOP FOR EACH WINDOW #################
	##############################################################
	
	#window sizes
	window_sizes = ['50kb', '100kb', '200kb', '500kb', '1000kb']
	
	#for each window size
	for i in range(0, len(window_sizes)):
		
		#select the [i] window
		selected_window = window_sizes[i]
			
		#open an empty list to save the dataframes subsetted
		list_dataframes = []
	
		#add each predictor path (these paths names were used to obtain both the data.frames and the files names)
		for k in range(0, len(list_predictors_files_full_path)):
	
			#select the [k] predictor name
			selected_predictor_name = list_predictors_files_only_names[k]
			
			#select the [k] data.frame
			selected_df = list_df_predictors[k]
			
			#list for save the variable names. We only include 'gene_id' for now.
			names_selected_predictor = ['gene_id']

			#if the selected predictor is NOT one of the predictors without windows
			if selected_predictor_name not in ['vip_distance_file.txt', 'pip_v3.txt', 'gene_length_v1.txt', 'gene_expression_final_v5.txt']:

				#extract the name of the column for the [k] window
				names_selected_predictor.append(''.join([selected_df.columns[1].split('50kb')[0], selected_window]))
			else: #if not
	
				#if the selected predictor names is NOT gene expression
				if selected_predictor_name not in ['gene_expression_final_v5.txt']:
	
					#extract the name of the column for the [k] window
					names_selected_predictor.append(selected_df.columns[1])
				else: #if not, and hence the selected predictor is gene expression
	
					#set mean gene expression as the predictor. For now, only:
					names_selected_predictor.append('mean_expres_all_tissues')
					names_selected_predictor.append('Testis')
					names_selected_predictor.append('Cells - EBV-transformed lymphocytes')
	
			#create a subset for gene id and the selected predictor
			subset_raw = selected_df[names_selected_predictor]

			#change the name of the second column using the file name (the first one is always gene_id)
			if any(check in names_selected_predictor for check in ['mean_expres_all_tissues', 'Testis', 'Cells - EBV-transformed lymphocytes']):
				
				#change the names of the columns by hand
				subset_raw = subset_raw.rename(columns = {'mean_expres_all_tissues': 'expression_all_tissues', 'Testis': 'testis_expression', 'Cells - EBV-transformed lymphocytes': 'lympho_expression'})
			else: #if not
				#change the name of the second column using the file name (the first one is always gene_id)
				subset_raw = subset_raw.rename(columns = {subset_raw.columns[1]: '_'.join([re.split(', |_|.txt|!',selected_predictor_name)[0], re.split(', |_|.txt|!',selected_predictor_name)[1]])})
			
			#save in the selected column in the list of dataframes
			list_dataframes.append(subset_raw)

		#merge all the dataframes within the list using 'gene_id'
		final_predictors = reduce(lambda x, y: pd.merge(x, y, on='gene_id', how='outer'), list_dataframes)

		#select the mean iHS and the number of iHS data points of the selected window
		ihs_selected_window = ihs_data[['gene_id', ''.join(['mean_ihs_', selected_window])]]
		ihs_n_points_selected_window = ihs_n_points[['gene_id', ''.join(['n_ihs_', selected_window])]]

		#change the name of the column with the number of iHS data points
		ihs_n_points_selected_window = ihs_n_points_selected_window.rename(columns = {ihs_n_points_selected_window.columns[1]: '_'.join([re.split(', |_|!',ihs_n_points_selected_window.columns[1])[0], re.split(', |_|!',ihs_n_points_selected_window.columns[1])[1]])})

		#bind iHS data and the number of data points for the selected window
		final_ihs_data = pd.merge(ihs_selected_window, ihs_n_points_selected_window, on='gene_id', how='outer')

		#bind with the predictors
		final_data = pd.merge(final_ihs_data, final_predictors, on='gene_id', how='outer')

		#remove the NAs
		final_data = final_data.dropna(axis=0, how='any')

		
		
		#########################################################
		################ PREPARE THE FINAL DATA #################
		#########################################################	
		
		#extract the column names as a header
		header = np.array(final_data.columns)
		
		#select all predictors and iHS but not the first name, the gene_id
		header = header[1:len(header)]
		
		#select the all the columns except the gene id and then transform to numpy array.
		data = final_data[header].values 
		
		# log transform response (iHS) and scale
		data[:, 0] = preprocessing.scale(np.log(data[:, 0])) 

		# standardize the whole data, including iHS and the predictors
		data = preprocessing.scale(data) 

		

		#################################################################
		##### DEFINE THE OBJECT CLASS THAT WILL PERFORM THE ANALYSES ####
		#################################################################
				
		#set the a new class of object called MixtureDensityModel
		class MixtureDensityModel(object):
		
			########################################
			####### SET THE INITIAL ACTIONS ########
			########################################	
			def __init__(self, data, response_index=0, m0=0., m1=0.5, sd0=1., sd1=1.): 
			
				#create the mask for subseting the columns
				mask = np.ones(data.shape[1], dtype=bool)
				
				#set the first element of mask as '0' instead of True, i.e., False.
				mask[response_index] = 0
				
				#from data, take all the rows of the columns indicated as true in mask (only the first is False)
				self.X = data[:, mask]
				
				#from data select the column of the response index
				self.Y = data[:, response_index]
		
				# add dummy variable (1) for bias term
				dummy_one = np.ones((data.shape[0], 1)) 

				#bind the column with the bias to the array of predictors
				self.X = np.concatenate((dummy_one, self.X), 1)
		
				#weights of covariates
				self.weight = np.zeros(self.X.shape[1])
		
				#parameters in the mixture model
				self.mixture_mean = np.array([m0, m1])
				self.mixture_sd = np.array([sd0, sd1])
				
				#the collection of parameters. This includes slopes of the predictors and Gaussian parameters (mean and sd).
				para = np.concatenate((self.mixture_mean, self.mixture_sd, self.weight)) 
		
				#set the object disable weight
				self.disabled_weight = None
		
				#initialize log likelihood
				self.eval(para)


		
			###################################
			####### SET THE EVALUATION ########
			###################################
			def eval(self, para):
		
				#save the parameters as an object called para inside the object class MixtureDensityModel
				self.para = para 
				
				#now do the same with the mean, sd and the weights extracting them from self.para
				self.mixture_mean = self.para[0:2]
		
				#the 2-3 elements of 'para' are the two sd.
				self.mixture_sd = self.para[2:4]
				
				#the 4 to the end elements of 'para' are the weights of the covariates including Bias (intercept).
				self.weight = self.para[4:]
		
				#multiply the values of each covariate by the weights
				xw = np.dot(self.X, self.weight)
					
				#save xw dot product
				self.expected_ihs_no_prob = xw

				#Sigmoid function to have values between 0 and 1.
				p = expit(xw)
 
				#calculate the density distribution of iHS for each gaussian distribution
				density0 = sp.stats.norm.pdf(self.Y, loc=self.mixture_mean[0], scale=self.mixture_sd[0])
				density1 = sp.stats.norm.pdf(self.Y, loc=self.mixture_mean[1], scale=self.mixture_sd[1])

				#now calculate the log likelihood
				self.log_likelihood = np.sum(np.log(density0 * (1 - p) + density1 * p)) 
		
				#save the point probability as p 
				self.point_prob = p 
			
		
		
			########################
			####### SET FIT ########
			########################
			def fit(self, likelihood_ratio_test=True):
				
				# estimate the best parameters (under the full model, .i.e., with all covariates) 
				self._fit(self.para)
		
				#perform the likelihood ratio test for each covariate
				if likelihood_ratio_test:
					self.p_value = self.likelihood_ratio_test()
			
		
		
			#########################
			####### SET _FIT ########
			#########################
			def _fit(self, para):
		
				#Using a list of parameters, we apply the function optimize.minimize to optimize 'para'
				res = optimize.minimize(self, para, method='L-BFGS-B') 
		
				#use the optimized parameters (slopes and Gaussian parameters) for running def __call__
				self(res.x) 
				
		
		
			#############################
			####### SET __call__ ########
			#############################
			def __call__(self, para):	

				#copy the parameters in all_para. This includes slopes of the covariates and parameters of the Gaussian distributions.
				all_para = np.copy(para)
		
				#if self.disabled_weight has a number, you remove the corresponding weight.
				if self.disabled_weight is not None:
					
					#insert a zero value in the position indicated by the index self.disabled_weight
					all_para = np.insert(all_para, self.disabled_weight, 0.) 
				
				#run eval function with all_para
				self.eval(all_para) 
		
				#return the - log_likelihood obtained from running the eval function
				return -self.log_likelihood 
			
		
		
			##########################################
			####### SET LIKELIHOOD_RATIO_TEST ########
			##########################################
			def likelihood_ratio_test(self):
				
				#copy the parameters stores in self.para, which came from the first optimization on the dummy parameters (def __init__), so we already have optimized parameters. All covariates were included. This will be called original parameters.
				ori_para = np.copy(self.para)
		
				#save the likelihood.
				ori_likelihood = self.log_likelihood
		
				#create an empty list to save the p.values 
				p_value_list = list()
					
				#loop for doing the likelihood ratio tests
				for para_idx in range(4, ori_para.shape[0]):

					#save the index of the covariate that will be disabled (weight remove) in this iteration
					self.disabled_weight = para_idx
					
					#copy the original parameters in a new object called para
					para = np.copy(ori_para)
		
					#from para, delete the weight of the covariate that correspond to the index (para_idx)
					para = np.delete(para, para_idx)
					
					#run def __call__ with the remaining parameters
					self(para) 
					
					#run self._fit to make the optimization of the parameters.
					self._fit(para)
		
					#likelihood ratio statistic
					lln_ratio = 2. * (ori_likelihood - self.log_likelihood) 
		
					# follow chi square distribution
					p_value = 1 - sp.stats.chi2.cdf(lln_ratio, 1) 
		
					#append the p vale into the previously empty list created
					p_value_list.append(p_value)
		
				#set disabled_weight to none
				self.disabled_weight = None 
				
				#This line is done to run def __call__ and hence run one last time the eval function, this time with the final means, sds, and weights, calculating in that way the probability of adaptation according the weighted input of ALL covariates and the corresponding distributions and mixture. 
				self(ori_para)
		
				#return the list of pvalues
				return np.array(p_value_list)
		
		
		
		##########################################
		########## RUN THE ANALYSES ##############
		##########################################
		
		#set the model using the MixtureDensityModel class applied to data
		model = MixtureDensityModel(data)
		
		#use the fit function that is defined inside the MixtureDensityModel class. Therefore, this function will be applied to the data included in model object
		model.fit(True)
		
		
		##save slopes and p-values
		#extract the header names.
		name = np.copy(header) #copy header names
		name[0] = "Intercept" #set the first name as Intercept.
		
		#create a panda data frame with the wieght and pvalues
		para_table = pd.DataFrame({'covariate': name, 'slope': model.weight, 'p_value': model.p_value})
		

		##change the names of the covariates	
		#create two lists with the original and extended covariate names
		raw_names = ["Intercept", "n_ihs", "vip_distance", "expression_all_tissues", "testis_expression", "lympho_expression", "chip_immune", "cons_elements", "pip_v3", "gene_length", "gene_number", "chip_testis", "tbfs_density", "coding_density", "gc_content", "recombination_final", "chip_density"]
		new_names = ["Intercept", "Number iHS data points", "Distance to VIPs", "Gene expression", "Gene expression in testis", "Gene expression in immune cells", "Regulatory density in immune cells (ChIP-seq)", "Density of conserved elements", "Number PPIs", "Gene length", "Gene number", "Regulatory density in testis (ChIP-seq)", "Regulatory density (DNaseI)", "Coding density", "GC-content", "Recombination rate", "Regulatory density (ChIP-seq)"]

		#replace the the values of covariates with the new names
		para_table.iloc[:,0] = para_table.iloc[:,0].replace(to_replace=raw_names, value=new_names)

		#save the table
		para_table.to_csv(path_or_buf=path_outside_container_outputs + '/tables' + '/' + pop_abbreviations[selected_pop] + '_' + selected_window + '_slopes_pvalues.csv', sep=',', header=True, index=False)

		
		##plot the distribution
		#open the plot
		axes = plt.gca() #this only plots the axes

		#set the font size for the whole plot
		plt.rcParams.update({'font.size': 13})
		
		#set the limits of the x axis
		axes.set_xlim([-5, 5])

		#plot the histogram of the response, which is iHS
		plt.hist(model.Y, 50, density=True, alpha=0.4, label='Observed iHS')
		
		#create distributions for the two components
		#define a function to do that
		def plot_gaussian(mean, sd, xmin, xmax, weight, label, shape):
			
			#make a range betwwen the lowest and highest value of x (iHS), being each step of 0,01
			x = np.arange(xmin, xmax, 0.01)

			#the y values would be a distribution 
			y = sp.stats.norm.pdf(x, mean, sd) * weight

			#plot using x, y, the shape of the line, label, width of the line and intensity of the color with alpha
			plt.plot(x, y, shape, label=label, linewidth=2, alpha=0.7)
		
			#return the values of x and y
			return x, y
		
		#extract the values of x (iHS) and y (density of the two components)
		x, y1 = plot_gaussian(mean=model.mixture_mean[0], sd=model.mixture_sd[0], xmin=-5., xmax=5., weight=1. - np.mean(model.point_prob), label='Component 1' , shape='r-') 
		x, y2 = plot_gaussian(mean=model.mixture_mean[1], sd=model.mixture_sd[1], xmin=-5., xmax=5., weight=np.mean(model.point_prob), label='Component 2', shape='r--') 
				
		#plot also the mixture of both components
		plt.plot(x, y1 + y2, 'b-', linewidth=3, label='Mixture', alpha=0.7) 
		
		#plot the xlabel
		plt.xlabel("log(iHS)")
		
		#plot the legend upper right
		plt.legend(loc='upper right')
		
		#plot the title
		plt.title('Population: ' + selected_pop + ' - Window size: ' + selected_window)
		
		#we use instead plt.savefig to save the figure
		plt.savefig(path_outside_container_outputs + '/figures' + '/' + pop_abbreviations[selected_pop] + '_' + selected_window + '_final_distributions.png', bbox_inches='tight')

		#close the plot
		plt.close()



############################################################
################## PARALLELIZE THE FUNCTION ################
############################################################

#open a pool
pool = mp.Pool(processes=5) #we set 5 processes.

#apply the function across the pool
results_pararllelize = [pool.apply_async(ihs_modelling_per_pop, args=(selected_pop,)) for selected_pop in list_pops_unique]

#close the pool
pool.close()

#wait for the completion of all scheduled jobs
pool.join()


##check we have all the expected files
#calculate the number of tables and figures
n_files = int(subprocess.check_output("find " + path_outside_container_outputs + " -type f | wc -l", shell=True))

#check
print('##############################################'); print('CHECK WE HAVE THE CORRECT NUMBER OF TABLES AND FIGURES AFTER RUNNING ANALYSES'); print('##############################################') 
print(n_files == 50)
print('##############################################')
