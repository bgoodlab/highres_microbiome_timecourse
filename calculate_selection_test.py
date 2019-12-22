import numpy
from numpy.random import poisson, binomial, beta, choice
from math import exp,log10,log
from scipy.special import gammaln
import sys
import pylab
from montecarlo_inference import *
from deterministic_inference import *
from deterministic_inference import infer_deterministic_trajectory, calculate_marginalized_overfit_loglikelihoods, infer_all_deterministic_intervals, calculate_number_nonfixed_intervals
from neutral_hmm import calculate_hmm_loglikelihoods, calculate_interval_log_pvalues, calculate_interval_zvalues
from neutral_hmm import *
        

if __name__=='__main__':
    
    data = []
    
    data = []
      
    
    # Alistipes example
    # Barcode version 07/19
    data_str = "0,27,474; 4,37,293; 10,0,39; 25,6,154; 35,28,250; 37,4,121; 46,1,269; 49,1,91; 67,305,379; 69,72,81; 72,299,392; 74,40,58; 77,80,132; 80,110,265; 92,30,85; 110,52,177; 118,80,139; 129,193,395; 153,114,262"
    # Barcode version w/ cluster without previous timepoint
    #data_str = "0,39,618; 4,44,404; 10,0,53; 25,7,191; 35,31,314; 37,6,169; 46,3,343; 49,1,117; 67,390,481; 69,94,108; 72,403,528; 74,44,68; 77,96,160; 80,159,374; 92,35,112; 110,71,229; 118,97,171; 129,247,501; 153,139,338"
    data.append(("Alistipes_onderdonkii", data_str, 0.75))
    
    # Phascolactobacterium example    
    data_str = "0,7,79; 4,292,508; 10,5,87; 25,25,255; 35,29,524; 37,11,154; 46,90,1324; 49,13,155; 67,485,485; 69,125,125; 72,548,548; 74,718,720; 77,197,197; 80,884,926; 92,258,270; 110,124,152; 118,74,182; 153,40,878"
    data.append(("Phascolactobacterium", data_str, 0.75))
    
    # Eligens
    # Barcode version 07/19
    data_str = "0,86404,913825; 4,195363,3184506; 10,25182,113935; 25,41245,1060154; 35,126238,1675794; 37,28295,341523; 46,55339,1430975; 49,18434,202148; 67,2031,1169754; 69,1638,290775; 72,96881,1342634; 74,1602022,1824237; 77,248345,354085; 80,5081350,5517747; 92,64222,126111; 110,133081,1268394; 118,110330,299425; 129,15033,242884; 153,48007,195101"
    data.append(("Eubacterium_eligens", data_str, 0.75))
    
    
    
    for species_name, data_str, starting_freq_threshold in data:
        print ""
        print species_name
        ts = []
        As = []
        Ds = []
        for item in data_str.split(";"):
            subitems = item.split(",")
            ts.append(long(subitems[0]))
            As.append(long(subitems[1]))
            Ds.append(long(subitems[2]))
    
        days = numpy.array(ts)
        As = numpy.array(As)
        Ds = numpy.array(Ds)
        ts = numpy.array(days,dtype=numpy.int32)
        fs = (As)*1.0/(Ds)

        clipped_Ds = numpy.clip(Ds,0,2000)
        clipped_As = numpy.array(numpy.around(fs*clipped_Ds),dtype=numpy.int32)
        
        As = clipped_As
        Ds = clipped_Ds
        fs = (As)*1.0/(Ds)
        
        print ts
        print As
        print Ds
        print fs
        
        birth_idx, last_low_idx, first_high_idx, death_idx = calculate_trajectory_milestones(ts,As,Ds)
        
        print "Birth:", birth_idx
        print "Death:", death_idx
        print "20-80:", last_low_idx, first_high_idx
    
        print "Jorde Ryman scores"
        Fprimes = calculate_jorde_ryman_scores(ts,As,Ds)
        print Fprimes
        print "All Ne:", 1.0/Fprimes.mean()
        
        Fprimes = calculate_jorde_ryman_scores(ts[first_high_idx:],As[first_high_idx:],Ds[first_high_idx:])
        print "Post-fixation Ne:", 1.0/Fprimes.mean()
        
        Fprimes = calculate_jorde_ryman_scores(ts[:last_low_idx+1],As[:last_low_idx+1],Ds[:last_low_idx+1])
        print "Pre-fixation Ne:", 1.0/Fprimes.mean()
        
        tfix = ts[first_high_idx]-ts[last_low_idx]
        print "Fixation time:", tfix
        
        taus,loglikelihoods = calculate_most_likely_time(As[last_low_idx]+1,Ds[last_low_idx]+2,As[first_high_idx]+1,Ds[first_high_idx]+2)
        print "Fixation Ne:", tfix/taus[loglikelihoods.argmax()]
        tauhat = taus[loglikelihoods.argmax()] 
        Nhat = 100.0
        #print "fs:", fs[last_low_idx],fs[first_high_idx]
        #print "As:", As[last_low_idx],As[first_high_idx]
        #print "Ds:", Ds[last_low_idx],Ds[first_high_idx]
        
        effective_growth_rate = numpy.ceil(tauhat*Nhat/tfix)
        Nhat = long(effective_growth_rate*tfix/tauhat)
        print "Effective generations per day:", effective_growth_rate
        print "Nhat:", Nhat
        print ts*effective_growth_rate
        ts = ts*effective_growth_rate
        

        null_Ns = numpy.logspace(log10(Nhat),log10(Nhat)+2,21)

        post_fixation_As = numpy.array(As)
        post_fixation_Ds = numpy.array(Ds)
        post_fixation_As[:first_high_idx] = 0
        post_fixation_Ds[:first_high_idx] = 0
        
        pre_fixation_As = numpy.array(As)
        pre_fixation_Ds = numpy.array(Ds)
        pre_fixation_As[last_low_idx+1:] = 0
        pre_fixation_Ds[last_low_idx+1:] = 0
        
        Ass = [As,post_fixation_As,pre_fixation_As]
        Dss = [Ds,post_fixation_Ds,pre_fixation_Ds]
        pivot_t_idxs = [first_high_idx, first_high_idx, last_low_idx]
        
        null_loglikelihoods = []
        for Ne in null_Ns:
            print "Ne =", Ne
            loglikelihoods = calculate_hmm_loglikelihoods(ts,Ass,Dss,Ne,pivot_t_idxs)
            null_loglikelihoods.append(loglikelihoods)
        
        null_loglikelihoods = numpy.array(null_loglikelihoods)
        
        all_Ne = null_Ns[null_loglikelihoods[:,0].argmax()]
        post_Ne = null_Ns[null_loglikelihoods[:,1].argmax()]
        pre_Ne = null_Ns[null_loglikelihoods[:,2].argmax()]
        postpre_Ne = null_Ns[(null_loglikelihoods[:,1]+null_loglikelihoods[:,2]).argmax()]
        
        
        print "All Ne:", all_Ne
        print "Post Ne:", post_Ne
        print "Pre Ne:", pre_Ne
        print "Postpre Ne:", postpre_Ne
    
        Nhat = postpre_Ne
        
        sys.stderr.write("Generating bootstrap samples!\n")
        num_bootstraps = 10000
        
        
        pre_f0s = sample_frequencies_from_data(ts,As,Ds,last_low_idx,size=num_bootstraps)
        post_f0s = sample_frequencies_from_data(ts,As,Ds,first_high_idx,size=num_bootstraps)
        
        post_ts = ts[first_high_idx:]-ts[first_high_idx]
        post_As = As[first_high_idx:]
        post_Ds = Ds[first_high_idx:]
        
        pre_ts = (ts[last_low_idx]-ts[:last_low_idx+1])[::-1]
        pre_As = (As[:last_low_idx+1])[::-1]
        pre_Ds = (Ds[:last_low_idx+1])[::-1]
        
        Ne_map = {'all': (all_Ne, all_Ne), 'joint': (postpre_Ne, postpre_Ne), 'separate': (post_Ne, pre_Ne)}   
        
        for key in sorted(Ne_map):
            
            print key
            
            Ne1, Ne2 = Ne_map[key]
        
            bootstrapped_post_fss = simulate_frequency_trajectories(post_ts,0,Ne1,post_f0s)
            bootstrapped_post_Ass = simulate_read_trajectories(post_ts, post_As, post_Ds,bootstrapped_post_fss)
            bootstrapped_pre_fss = simulate_frequency_trajectories(pre_ts,0,Ne2,pre_f0s)
            bootstrapped_pre_Ass = simulate_read_trajectories(pre_ts, pre_As, pre_Ds,bootstrapped_pre_fss)
        
            sys.stderr.write("Done!\n")
        
            sys.stderr.write("Calculating post HMM likelihoods...\n")
            Ass = numpy.vstack([post_As, bootstrapped_post_Ass])
            Dss = numpy.vstack([post_Ds, numpy.array([post_Ds]*num_bootstraps)])
        
            loglikelihoods = calculate_hmm_loglikelihoods(post_ts,Ass,Dss,Ne1)
        
            observed_post_loglikelihood = loglikelihoods[0]
            bootstrapped_post_loglikelihoods = loglikelihoods[1:]
            sys.stderr.write("Done!\n")
        
            observed_LRT = -observed_post_loglikelihood
            bootstrapped_LRTs = -bootstrapped_post_loglikelihoods
            pvalue = ((bootstrapped_LRTs>=observed_LRT).sum()+1.0)/(len(bootstrapped_LRTs)+1.0)
        
            print "Post pvalue:", pvalue
        
            sys.stderr.write("Calculating pre HMM likelihoods...\n")
            Ass = numpy.vstack([pre_As, bootstrapped_pre_Ass])
            Dss = numpy.vstack([pre_Ds, numpy.array([pre_Ds]*num_bootstraps)])
        
            loglikelihoods = calculate_hmm_loglikelihoods(pre_ts,Ass,Dss,Ne2)
        
            observed_pre_loglikelihood = loglikelihoods[0]
            bootstrapped_pre_loglikelihoods = loglikelihoods[1:]
        
            sys.stderr.write("Done!\n")
        
            observed_LRT = -observed_pre_loglikelihood
            bootstrapped_LRTs = -bootstrapped_pre_loglikelihoods
            pvalue = ((bootstrapped_LRTs>=observed_LRT).sum()+1.0)/(len(bootstrapped_LRTs)+1.0)
        
            print "Pre pvalue:", pvalue
        
            observed_loglikelihood = observed_post_loglikelihood + observed_pre_loglikelihood
            bootstrapped_loglikelihoods = bootstrapped_pre_loglikelihoods + bootstrapped_post_loglikelihoods
        
            observed_LRT = -observed_loglikelihood
            bootstrapped_LRTs = -bootstrapped_loglikelihoods
            pvalue = ((bootstrapped_LRTs>=observed_LRT).sum()+1.0)/(len(bootstrapped_LRTs)+1.0)
        
            print "Both pvalue:", pvalue
        
    sys.exit(0)
    generation_time = 10
    
    for species_name, data_str, starting_freq_threshold in data:
    
        print species_name
        ts = []
        As = []
        Ds = []
        for item in data_str.split(";"):
            subitems = item.split(",")
            ts.append(long(subitems[0]))
            As.append(long(subitems[1]))
            Ds.append(long(subitems[2]))
    
        days = numpy.array(ts)
        As = numpy.array(As)
        Ds = numpy.array(Ds)
        ts = numpy.array(days,dtype=numpy.int32)*generation_time
        fs = (As)*1.0/(Ds)

        clipped_Ds = numpy.clip(Ds,0,2000)
        clipped_As = numpy.array(numpy.around(fs*clipped_Ds),dtype=numpy.int32)
        
        As = clipped_As
        Ds = clipped_Ds
        fs = (As)*1.0/(Ds)
        
        birth_idx, last_low_idx, first_high_idx, death_idx = calculate_trajectory_milestones(ts,As,Ds)
        
        print birth_idx, last_low_idx, first_high_idx, death_idx
    
        print "Jorde Ryman scores"
        Fprimes = calculate_jorde_ryman_scores(ts,As,Ds)
        print Fprimes
        print 1.0/Fprimes.mean()
        
        print "Post-fixation"
        Fprimes = calculate_jorde_ryman_scores(ts[first_high_idx:],As[first_high_idx:],Ds[first_high_idx:])
        print 1.0/Fprimes.mean()
        
        print "Pre-fixation"
        Fprimes = calculate_jorde_ryman_scores(ts[:last_low_idx+1],As[:last_low_idx+1],Ds[:last_low_idx+1])
        print 1.0/Fprimes.mean()
        
        print "Fixation"
        tfix = ts[first_high_idx]-ts[last_low_idx]
        print tfix
        
        print "Calculating better value..."
        taus,loglikelihoods = calculate_most_likely_time(As[last_low_idx],Ds[last_low_idx],As[first_high_idx],Ds[first_high_idx])
        print taus[loglikelihoods.argmax()]
        print tfix/taus[loglikelihoods.argmax()]
        print "Done!"
        
        print taus
        print loglikelihoods
        
        
        if starting_freq_threshold>-0.5:
            #start_idx = 0
            
            #start_idx = numpy.nonzero(ts>60*generation_time)[0][0]
            
            start_idx = numpy.nonzero(fs>=starting_freq_threshold)[0][0]
            
            ts = ts[start_idx:]
            As = As[start_idx:]
            Ds = Ds[start_idx:]
            fs = fs[start_idx:]
    
            pivot_t_idx = 0
        else:
            
            pivot_t_idx = (fs*(1-fs)).argmax()
    
        post_fixation_ts = ts[first_high_idx:]
        post_fixation_As = As[first_high_idx:]
        post_fixation_Ds = Ds[first_high_idx:]
        post_fixation_pivot_t_idx = 0
        
        
    
        print pivot_t_idx
        print ts
        print As
        print Ds
        print fs
    
        
        Ys = calculate_fit_scores(ts,As,Ds)
        print "Fit scores"
        print Ys
        print Ys.mean(), 1.0/numpy.square(Ys).mean()
        print Ys.std()
    
        n = calculate_number_nonfixed_intervals(ts,As,Ds)
    
        sys.stderr.write("Running logistic regression...\n")
        deterministic_s, deterministic_s_CI, f0, deterministic_loglikelihood = infer_deterministic_trajectory(ts,As,Ds)
        s_lower = deterministic_s-(deterministic_s-deterministic_s_CI[0])
        s_upper = deterministic_s+(deterministic_s_CI[1]-deterministic_s)
        print deterministic_s, deterministic_s_CI
        sys.stderr.write("Done!\n")
        
        
        sys.stderr.write("Calculating all logistic regression...\n")
        interval_s_map = infer_all_deterministic_intervals(ts,As,Ds)
        sys.stderr.write("Done!\n")
        
        #print interval_s_map
        
        sys.stderr.write("Pre-computing null trajectories...\n")
        #null_Ns = numpy.logspace(2,4,3)
        null_Ns = numpy.logspace(2,4,11)
        #null_Ns = numpy.array([100,1000,10000])
        null_trajectories = []
        null_trajectory_loglikelihoods = []
        null_loglikelihoods = []
        for N in null_Ns:
            s=0
            sys.stderr.write("N=%g\n" % N)
            fss, loglikelihoods = simulate_model(ts,As,Ds,s,N,size=100000,starting_t_idx=pivot_t_idx)
            null_trajectories.append( fss )
            null_trajectory_loglikelihoods.append( loglikelihoods )
            null_loglikelihoods.append( calculate_average_loglikelihood(loglikelihoods) )
            
                
        null_trajectories = numpy.array(null_trajectories)
        null_trajectory_loglikelihoods = numpy.array(null_trajectory_loglikelihoods)
        null_loglikelihoods = numpy.array(null_loglikelihoods)
        sys.stderr.write("Done!\n")
    
        
        alternative_params = []
        #alternative_Ns = numpy.logspace(3,6,4)
        alternative_Ns = [1e05]
        for N in alternative_Ns:
            alternative_params.append((deterministic_s,N))
        #ss = numpy.logspace(log10(s_lower),log10(s_upper),5)
        #ss = numpy.hstack([[deterministic_s],ss])
        #for s in ss:
        #    N = 1e06 # Deterministic
        #    alternative_params.append((s,N))
    
        sys.stderr.write("Pre-computing alternative trajectories...\n")
        alternative_trajectories = []
        alternative_trajectory_loglikelihoods = []
        alternative_loglikelihoods = []
        for s,N in alternative_params:
            sys.stderr.write("s=%g, N=%g\n" % (s,N))
            
            fss, loglikelihoods = simulate_model(ts,As,Ds,s,N,size=100000,starting_t_idx=pivot_t_idx)
            
            alternative_trajectories.append( fss )
            alternative_trajectory_loglikelihoods.append( loglikelihoods ) 
            alternative_loglikelihoods.append( calculate_average_loglikelihood(loglikelihoods) )
            
        alternative_trajectories = numpy.array(alternative_trajectories)
        alternative_trajectory_loglikelihoods = numpy.array(alternative_trajectory_loglikelihoods)
        alternative_loglikelihoods = numpy.array(alternative_loglikelihoods)
        sys.stderr.write("Done!\n")
    
        # Find maximum likelihood estimates of neutral model
        print "Null loglikelihoods", null_loglikelihoods
        max_null_L = null_loglikelihoods.max()
        max_null_idx = null_loglikelihoods.argmax()
        
        # Find maximum likelihood estimates of deterministic selection model
        max_alt_L = alternative_loglikelihoods.max()
        max_alt_idx = alternative_loglikelihoods.argmax()
        
        print "Alternative loglikelihoods", alternative_loglikelihoods
        
        print "Null hypothesis", max_null_L, null_Ns[max_null_idx]
        print "Alternative hypothesis", max_alt_L, alternative_params[max_alt_idx]
        print "Deterministic alternative", deterministic_loglikelihood, deterministic_s
        observed_loglikelihood_ratio = max([0,max_alt_L-max_null_L])
        observed_L = max_null_L
    
        print "Observed loglikelihood ratio =", observed_loglikelihood_ratio
    
        Nhat = null_Ns[max_null_idx]
        s,Nhat
    
        sys.stderr.write("Bootstrapping...\n")
        num_bootstraps = 10000
        
        
        sys.stderr.write("Simulating bootstrapped trajectories...\n")
        posterior_fss = sample_posterior_trajectories(null_trajectories[max_null_idx],null_trajectory_loglikelihoods[max_null_idx],size=num_bootstraps)
        
        # Pure posterior distribution (not quite right...)
        # bootstrapped_Ass = simulate_read_trajectories(ts,As,Ds,posterior_fss)
        
        # Doing it around a pivot
        #pivot_t_idx = numpy.nonzero((fs>0.1))[0][0]
        #pivot_t_idx = (fs*(1-fs)).argmax()
        f0s = posterior_fss[:,pivot_t_idx]
        
        f0s = sample_frequencies_from_data(ts,As,Ds,pivot_t_idx,size=num_bootstraps)
        
        bootstrapped_fss = simulate_frequency_trajectories(ts,0,Nhat,f0s,pivot_t_idx)
        bootstrapped_Ass = simulate_read_trajectories(ts,As,Ds,bootstrapped_fss,pivot_t_idx)
        
        good_bootstrap_idxs = ((bootstrapped_Ass>0)*(bootstrapped_Ass<Ds[None,:])).any(axis=1)
        bootstrapped_Ass = bootstrapped_Ass[good_bootstrap_idxs,:]
        
        sys.stderr.write("Done!\n")
        
        sys.stderr.write("Calculating deterministic likelihoods...\n")
        d1,d2,d3,observed_deterministic_loglikelihood = infer_deterministic_trajectory(ts,As,Ds)
        observed_overfit_loglikelihood = calculate_marginalized_overfit_loglikelihoods(ts,As,Ds)
        observed_n = calculate_number_nonfixed_intervals(ts,As,Ds)
        
        bootstrapped_deterministic_loglikelihoods = []
        bootstrapped_overfit_loglikelihoods = []
        bootstrapped_ns = []
        for idx in xrange(0,bootstrapped_Ass.shape[0]):
        
            d1,d2,d3,deterministic_loglikelihood = infer_deterministic_trajectory(ts,bootstrapped_Ass[idx,:],Ds)
            bootstrapped_deterministic_loglikelihoods.append(deterministic_loglikelihood)
            
            overfit_loglikelihood = calculate_marginalized_overfit_loglikelihoods(ts,bootstrapped_Ass[idx,:],Ds)
            bootstrapped_overfit_loglikelihoods.append(overfit_loglikelihood)
            
            n = calculate_number_nonfixed_intervals(ts,bootstrapped_Ass[idx,:],Ds)
            bootstrapped_ns.append(n)
            
        bootstrapped_deterministic_loglikelihoods = numpy.array(bootstrapped_deterministic_loglikelihoods)    
        bootstrapped_overfit_loglikelihoods = numpy.array(bootstrapped_overfit_loglikelihoods)
        bootstrapped_ns = numpy.array(bootstrapped_ns)
        sys.stderr.write("Done!\n")
        
        #sys.stderr.write("Calculating interval pvalues...\n")
        #observed_log_ps, bootstrapped_log_ps = calculate_interval_log_pvalues(ts,As,Ds,Nhat,bootstrapped_Ass)
        #sys.stderr.write("Done!\n")
        
        #print observed_log_ps
        #print bootstrapped_log_ps[0:10]
        
        #observed_mean = -1*observed_log_ps.mean()
        #bootstrapped_means = []
        #for i in xrange(0,len(bootstrapped_log_ps)):
        #    bootstrapped_means.append(-1*bootstrapped_log_ps[i].mean())
        #bootstrapped_means = numpy.array(bootstrapped_means)
        
        #interval_pvalue = ((bootstrapped_means>=observed_mean).sum()+1.0)/(len(bootstrapped_means)+1.0)
        
        #print "interval pvalue =", interval_pvalue       
        
        sys.stderr.write("Calculating interval zscores...\n")
        observed_z_vector, bootstrapped_z_vector = calculate_interval_zvalues(ts,As,Ds,Nhat,bootstrapped_Ass)
        sys.stderr.write("Done!\n")
        
        print observed_z_vector
        print bootstrapped_z_vector[0:10]
        
        sys.stderr.write("Calculating HMM likelihoods...\n")
        observed_loglikelihood, bootstrapped_loglikelihoods = calculate_hmm_loglikelihoods(ts,As,Ds,Nhat,bootstrapped_Ass)
        sys.stderr.write("Done!\n")
        
        hmm_pvalue = (((bootstrapped_loglikelihoods/bootstrapped_ns)<=(observed_loglikelihood/observed_n)).sum()+1.0)/(len(bootstrapped_loglikelihoods)+1.0)

        print "HMM loglikelihood =", observed_loglikelihood
        print "HMM pvalue =", hmm_pvalue
        
        
        observed_LRT = max([0,observed_deterministic_loglikelihood-observed_loglikelihood])
        bootstrapped_LRTs = numpy.clip(bootstrapped_deterministic_loglikelihoods-bootstrapped_loglikelihoods,0,1e09)
        
        LRT_pvalue = ((bootstrapped_LRTs>=observed_LRT).sum()+1.0)/(len(bootstrapped_LRTs)+1.0)
        
        print "LRT pvalue =", LRT_pvalue
        
        observed_LRT = (observed_overfit_loglikelihood-observed_loglikelihood)/observed_n
        bootstrapped_LRTs = (bootstrapped_overfit_loglikelihoods-bootstrapped_loglikelihoods)/bootstrapped_ns
        overfit_LRT_pvalue = ((bootstrapped_LRTs>=observed_LRT).sum()+1.0)/(len(bootstrapped_LRTs)+1.0)
        
        print "Overfit loglikelihood =", observed_overfit_loglikelihood
        print "Overfit LRT pvalue =", overfit_LRT_pvalue
        
        sys.stderr.write("Plotting comparison...\n")
        
        # Plot fitted trajectories
        pylab.figure(figsize=(7,9))
        #pylab.xlim([0,1010])
        
        pylab.subplot(311)
        pylab.title("%s, Plrt=%g, Pnew=%g" % (species_name, LRT_pvalue, hmm_pvalue))
        
        pylab.xlim([ts[0]/Nhat,ts[-1]/Nhat])
        pylab.ylim([0,1])
        pylab.plot(ts/Nhat,fs,'ko')
        pylab.ylabel('Allele freq')
        
        for i in xrange(0,20):
            pylab.plot(ts/Nhat,bootstrapped_fss[i,:],'-')
        
        
        pylab.subplot(312)
        pylab.xlim([ts[0]/Nhat,ts[-1]/Nhat])
        pylab.ylim([0,1])
        pylab.plot(ts/Nhat,fs,'ko')
        pylab.ylabel('Allele freq')
        
        
        posterior_fss = sample_posterior_trajectories(null_trajectories[max_null_idx],null_trajectory_loglikelihoods[max_null_idx],size=20)
        for i in xrange(0,posterior_fss.shape[0]):
            pylab.plot(ts/Nhat,posterior_fss[i],'-')
        
        pylab.subplot(313)
        pylab.xlim([ts[0]/Nhat,ts[-1]/Nhat])
        pylab.ylim([0,1])
        pylab.plot(ts/Nhat,fs,'ko')
        pylab.xlabel('Time (generations)')
        pylab.ylabel('Allele freq')
        posterior_fss = sample_posterior_trajectories(alternative_trajectories[max_alt_idx],alternative_trajectory_loglikelihoods[max_alt_idx],size=20)
        for i in xrange(0,posterior_fss.shape[0]):
            pylab.plot(ts/Nhat,posterior_fss[i],'-')
        filename = "%s.pdf" % (species_name)
        pylab.savefig(filename,bbox_inches='tight')
        sys.stderr.write("Done!\n")
        
        
        # How to do it numerically!    
        #bootstrapped_loglikelihood_ratios = []
        #bootstrapped_Ls = []
        #for fs in null_trajectories[max_null_idx][0:10]:
        
            # Now evaluate loglikelihoods
            #bootstrapped_null_loglikelihoods = []
            #for fss in null_trajectories:
            #    bootstrapped_null_loglikelihoods.append( evaluate_loglikelihood(ts,bootstrapped_As,Ds,fss) )
            #bootstrapped_null_loglikelihoods = numpy.array(bootstrapped_null_loglikelihoods)
    
            #bootstrapped_max_null_L = bootstrapped_null_loglikelihoods.max()
        
            #bootstrapped_Ls.append( bootstrapped_max_null_L )
        
            # Now evaluate loglikelihoods
            #bootstrapped_alternative_loglikelihoods = []
            #for fss in alternative_trajectories:
            #    bootstrapped_alternative_loglikelihoods.append( evaluate_loglikelihood(ts,bootstrapped_As,Ds,fss) )
            #bootstrapped_alternative_loglikelihoods = numpy.array(bootstrapped_alternative_loglikelihoods)
    
            #bootstrapped_max_alt_L = bootstrapped_alternative_loglikelihoods.max()
            #bootstrapped_loglikelihood_ratio = max([0,bootstrapped_max_alt_L-bootstrapped_max_null_L])
    
            #bootstrapped_loglikelihood_ratios.append( bootstrapped_loglikelihood_ratio )
            #print len(bootstrapped_loglikelihood_ratios), bootstrapped_loglikelihood_ratio, bootstrapped_max_null_L
        
            #if len(bootstrapped_Ls) % 100 == 0:
            #    sys.stderr.write("%d\n" % len(bootstrapped_Ls))
        
        #bootstrapped_loglikelihood_ratios = numpy.array(bootstrapped_loglikelihood_ratios)
        #bootstrapped_Ls = numpy.array(bootstrapped_Ls)
    
        #pvalue = ((bootstrapped_loglikelihood_ratios>=observed_loglikelihood_ratio).sum()+1.0)/(len(bootstrapped_loglikelihood_ratios)+1.0)   
        #print pvalue
    
        #pvalue = ((bootstrapped_Ls<=observed_L).sum()+1.0)/(len(bootstrapped_Ls)+1.0)   
        #print "Numerical pvalue =", pvalue
        
        sys.stdout.flush()
        
        