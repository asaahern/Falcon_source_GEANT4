#!/usr/bin/env python
import os
import time
import random
import multiprocessing  # the module we will be using for multiprocessing

def work(Run):

	WORKING_DIRECTORY=os.getcwd()
	os.system("./FALCON "+str(Run)+" 0")
	print ("Unit of work number %d" % Run )

if __name__ == "__main__":  # Allows for the safe importing of the main module
	print("There are %d CPUs on this machine" % multiprocessing.cpu_count())
	number_processes = multiprocessing.cpu_count()-1
	pool = multiprocessing.Pool(number_processes)
	total_tasks = 10
	tasks = range(total_tasks)
	results = pool.map_async(work, tasks)
	pool.close()
	pool.join()
	
