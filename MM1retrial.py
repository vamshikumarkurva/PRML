import numpy as np
import random

if __name__=="__main__":
	clock =0
	N = int(raw_input('enter the number of Telephone lines: '))
	AT = int(raw_input('enter mean time between two succesive arrivals(Arrival time): '))
	ST = int(raw_input('enter mean length of call duration(service time): '))
	RT = int(raw_input('enter the retrial time: '))
	LINKS = int(raw_input('enter the number of servers(Links): '))
	end_time = int(raw_input('enter the time for simulation: '))
	processed = 0
	busy = 0
	blocked = 0
	completed = 0
	links_in_use = 0
	status_lines = np.zeros(N)
	status_links = np.zeros(LINKS)
	lines_links = np.zeros(N)
	CAT = round(np.random.exponential(AT))  
	service_time = np.zeros(LINKS)
	retrial = np.zeros(N)
	retrial.fill(np.inf)
	waittime = np.zeros(N)
	SCT = np.zeros(N)
	SCT.fill(np.inf)

	while (clock < end_time):
	    not_busy_lines = np.where(status_lines==0)[0].tolist()
	    not_busy_links = np.where(status_links==0)[0].tolist()
	    min_time = min(SCT)
	    if ~all(retrial==np.inf):
		MIN = min(retrial);
		clock1 = min(CAT,min_time)
		clock1 = min(clock1,MIN)
	    else:
		clock1 = min(CAT,min_time)
	    
	    if (clock1==CAT or clock1==MIN):
	       if clock1==CAT:
		   CAT = CAT + round(np.random.exponential(AT))+1;
		   if len(not_busy_lines)>2:
		       call = np.random.choice(not_busy_lines,2,replace=False)
		   else:
		       clock = clock1
		       continue
	       elif clock1==MIN:
		   call = np.where(retrial==MIN)[0].tolist()
		   service_time = round(np.random.exponential(ST))
	       if (status_lines[call[0]]==1 or status_lines[call[1]]==1):
		   busy = busy+1
		   processed = processed+1
		   retrial[call] = clock1+RT
		   waittime[call] = waittime[call]+RT
		   status_lines[call] = 2
	       elif all(status_links==1):
		   blocked = blocked+1
		   processed = processed+1
		   retrial[call] = clock1+RT
		   waittime[call] = waittime[call]+RT
		   status_lines[call] = 2
	       else:
		   lines_links[call] = np.random.choice(not_busy_links,1)
		   status_lines[call] = 1 
		   status_links[lines_links[call[0]]]=1
		   SCT[call] = clock1+service_time
		   retrial[call] = np.inf
	    elif clock1==min_time:
		call = np.where(SCT==min_time)[0].tolist()
		status_lines[call] = 0
		status_links[lines_links[call[1]]]=0
		lines_links[call] = 0
		completed = completed + 1
		processed = processed + 1
		SCT[call] = np.inf 
	    clock = clock1
	
	avg_wait_time = sum(waittime)/(2*processed) 
	print 'Processed calls: ', processed
	print 'Blocked calls: ',blocked
	print 'Busy calls: ',busy
	print 'completed calls: ', completed
	print 'average waiting time',avg_wait_time
	
