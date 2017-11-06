

##### import functions
import pyemma
import pyemma.coordinates as coor
import pyemma.msm as msm
import pyemma.plots as mplt
from pyemma import config
import numpy as np
#import matplotlib.pyplot as plt
import pylab
import scipy
#import sklearn.cluster
import os
import subprocess as sbp
#import networkx as nx
#import pygraphviz as PG
import random as rd
import sys
from numpy.random import randint, permutation
from scipy.linalg import eigvals, inv
from pmx import *



class Test_the_model:
    
       
    '''
    Use this class, if you know the right structure and if you want tho evaluate your model quality.
    '''
    
    def __init__(self,reference, check_list):
        self.reference = reference
        self.check_list=check_list
        self.description = "Insert the refernece and the structures you want to compare. Possible formats are .pdb and .xtc"
        self.rmsd_list=[]
        
    def compare_to_reference(self, reference, out_data):
        self.rmsd_list=[]
        for entry in self.check_list:
            self.rmsd_list.append(Gromacs_commands(5).calc_rmsd_ref(reference, entry, out_data))
        return self.rmsd_list
    
    def rmsd_plot(self, class_structure, prob_score):
        pylab.plot( class_structure.prob[:5]*100, self.rmsd_list, 'o' )
        pylab.legend(bbox_to_anchor=(2.,1))
        pylab.grid(lw=2)
        if prob_score=='prob':
            pylab.xlabel('probability [%]', fontsize=16)
        elif prob_score=='score':
            pylab.xlabel('score', fontsize=16)
        pylab.ylabel('rmsd to reference [nm]', fontsize=16)
        pylab.ylim(0,3)
        pylab.show()
            


class Gromacs_commands:
    '''
    This class containes all commands for gromacs the program needs. Make sure to use the right gromacs version.
    '''
    def __init__(self,version):
        #gromacs version, insert eg. 4 for 4.6.7 or 5 for 5.2.3 !
        self.version = version
        self.description = "Please insert the gromacs version as eg. 4 for 4.6.7 or 5 for 5.2.3"
        
    def align_traj(self, ref, traj, out, sim_timestep ):
        if self.version==4:
                subs = sbp.Popen(['trj_conv','-s',ref, '-f', traj, '-dt', sim_timestep, '-o', out, '-fit', 'rot+trans' ], shell=False, stdin = sbp.PIPE)
        elif self.version==5:
                subs = sbp.Popen(['gmx', 'trjconv','-s',ref, '-f', traj, '-dt', sim_timestep, '-o', out, '-fit', 'rot+trans' ], shell=False, stdin = sbp.PIPE)
        print >> subs.stdin,"Protein"    # align proteins 
        print >> subs.stdin,"System"     # return whole system
        out,err = subs.communicate()
        
    def calc_mean_rmsd(self, outfile):
        f=open(outfile,'r')
        content=f.readlines()
        n=0
        rmsd_list=[]

        for n in range(len(content)):
            if content[n].startswith('#') or content[n].startswith('@'):
                pass
            else:
                rmsd_list.append(content[n])
        b=[]
        for m in range(len(rmsd_list)):
            a=str(rmsd_list[m]).split()
            b.append(a[1])
        c=np.array([b], dtype=float)
        rmsd_mean=np.mean(c)
        return rmsd_mean

        
    def calc_rmsd_ref(self, reference,  entry, out_data):
        outfile=out_data+'/rmsd.xvg'
        if self.version==4:
            subs = sbp.Popen(['g_rms','-s', reference, '-f', entry, '-o', outfile, '-fit', 'rot+trans' ], shell=False, stdin = sbp.PIPE)
        if self.version==5:
            subs = sbp.Popen(['gmx', 'rms','-s', reference, '-f', entry, '-o', outfile, '-fit', 'rot+trans' ], shell=False, stdin = sbp.PIPE)
        print >> subs.stdin,"C-alpha"
        print >> subs.stdin,"MOL"
        out,err = subs.communicate()
	try:
            rmsd=self.calc_mean_rmsd(outfile)
	except:
	    print 'rmsd calculation failed'
	    rmsd=1
        subs = sbp.Popen(['rm', out_data+'/rmsd.xvg'])
        return rmsd
    
    def calc_rmsd_traj(self, reference, traj_list, out_data):
        i=0
        rmsd_traj_list=[]
        for name in traj_list:
            outfile=out_data+'/rmsd%i.xvg' % (i)
            if self.version==4:
                subs = sbp.Popen(['g_rms','-s', reference, '-f',name, '-o', outfile, '-fit', 'rot+trans' ], shell=False, stdin = sbp.PIPE)
            if self.version==5:
                subs = sbp.Popen(['gmx', 'rms','-s', reference, '-f', name, '-o', outfile, '-fit', 'rot+trans' ], shell=False, stdin = sbp.PIPE)
            print >> subs.stdin,"C-alpha"
            print >> subs.stdin,"MOL"
            out,err = subs.communicate()
            i+=1
            rmsd_traj_list.append(outfile)
        return rmsd_traj_list
    
    def calc_b_factor(self, traj, input_pdb, rank, begin_time, end_time, out_data):
        outfile='%s/B_fac_%i.pdb' % (out_data, rank)
        if self.version==4:
            subs = sbp.Popen(['g_rmsf', '-f', traj,  '-s', input_pdb, '-nofit' , '-oq', outfile, '-b',  str(begin_time), '-e', str(end_time) ], shell=False, stdin = sbp.PIPE)
        if self.version==5:
            command_test='gmx rmsf -f %s -s  %s -nofit -oq %s/B_fac_%i.pdb -b %i -e %i' % (traj, input_pdb, out_data, rank, begin_time, end_time)
            subs = sbp.Popen(['gmx', 'rmsf', '-f', traj,  '-s', input_pdb, '-nofit' , '-oq', outfile, '-b',  str(begin_time), '-e', str(end_time) ], shell=False, stdin = sbp.PIPE)
            
        print >> subs.stdin,"MOL"
        out,err= subs.communicate()


class Load_process_data:
    '''
     This class loads the structure (.pdb file) and the simulations (.xtc files ) to build the Markov-state-model (MSM).
     Make sure that all .xtc files lay in the same folder. Also no additional .xtc files that do not belong to
     this MSM should be in there. For the .pdb file provide the whole path, for the .xtc files the folder in which all the
     .xtc file lay. The program will then search for files that end on .xtc
    '''
    def __init__(self, input_xtc, input_pdb, gromacs_version=5, sim_timestep=1, is_aligned=False, bootstrap=False):
        self.input_xtc = input_xtc
        self.input_pdb = input_pdb
        # returns the simulation timestep in nanoseconds
        self.sim_timestep= sim_timestep
        #states wether all simulations are aligned to the structure. Very important for heteroatom clustering 
        self.is_aligned= is_aligned
        # gives back, if user has chosen to do an bootstrap trial
        self.bootstrap= bootstrap
        self.traj_list = []
        self.gromacs_version=gromacs_version
        self.bootstrap_list=[]
        self.description = "This class loads the input structure and the trajectories in the desired form."
    
    def load_data(self):
        #loads the raw data
	for name in self.input_xtc:
        	try:
            		for filename in os.listdir(name):
                		if filename.endswith('.xtc'):
                    			self.traj_list.append(os.path.join(name,filename))
                    
        	except TypeError:
            		print "TypeError: Please enter the path as a string"
            
        	if not os.path.exists(name):
            		print 'InputError: Trajectory file does not excist'
            
        	if self.traj_list == []:
            		print "No .xtc trajectories were detected in your input folder"
            
        return self.traj_list
    
    def load_structure(self):
        self.topfile =  self.input_pdb
        #check if file excitst    
        if not os.path.isfile(self.topfile):
            print 'InputError: Structure file does not excist'
            
        return self.topfile
    
    def align_traj(self):
        #load first the structure, and the trajectory list
        self.topfile=self.load_structure()
        self.traj_list=self.load_data()
        #new trajectories will be saved in a new list
        self.outfile_list=[]
        i=0
        for name in self.traj_list:
            i+=1
            outfile=out_fit+'/fit%i.xtc' % i
            do_command=Gromacs_commands(self.gromacs_version)
            do_command.align_traj(self.topfile, name, outfile, str(self.sim_timestep))
            self.outfile_list.append(outfile) 
        return self.outfile_list
    
    def bootstrap_data(self):
        #first load the raw data
	self.bootstrap_list=[]
        #self.traj_list=self.load_data()
        len_list=len(self.traj_list)
        #randomly choose trajectories to sample (drawing with replacement)
        for k in range(len_list):
            rand_int=np.random.randint(len_list)
            self.bootstrap_list.append(self.traj_list[rand_int])
        return self.bootstrap_list


class Cluster_data:
    
    '''
    In this class first the coordinates on which to cluster have to be defined. You can cluster on the heteroatoms of the
    fragment (heteroatom), the smallest distance between the fragment`s and the protein`s atoms (distance) or on the 
    heteroatoms and important protein coordinate (prot_frag). The clustering methods are k_means clustering and 
    uniform_time clustering. 
    '''
    
    def __init__(self, traj_list, input_pdb, n_clusters, lagtime, cluster_method='uniform_time', coordinates='heteroatom', fragment_name='MOL', visual = False ):
        self.traj_list = traj_list
        self.input_pdb = input_pdb
        self.reduce_dimensions=False
        self.visual=visual
        self.dtrajs=0
        self.clustering_centers=0
        self.inp=0
        self.lagtime=lagtime
        self.n_clusters=n_clusters
        self.Y=0
        # Here the clustering method is defined
        self.cluster_method=  cluster_method
        self.coordinates= coordinates
        #fragment is per default thought to be named MOL
        self.fragment_name=fragment_name
        self.description = "Possible clustering methods are: k_means, uniform_time. Possible coordinates are: heteroatom, distance, prot_frag"
    
    def fetch_atoms(self):
        m = Model( self.input_pdb )
        res_names = self.fragment_name
        atom_names = ['N', 'N1', 'N2', 'N3', 'N4', 'NA', 'NB','ND', 'S', 'SS', 'S2', 'S1', 'O2','O', 'OS', 'O1','P','P1', 'P2', 'F', 'CL']  # list of possible heteroatoms than can occur in the fragments
        atom_names_backup = ['C', 'C1','C4','C5','C6', 'CA', 'CB', 'CC', 'CD', 'CE']
    
        # get atom indices
        res_list = m.fetch_residues( res_names )
        self.atom_ids = []
    
        for r in res_list:
            atoms = r.fetch_atoms( atom_names )
            for atom in atoms: self.atom_ids.append( atom.id - 1 )
            
        # check if at least 4 atomes were chosen if not, add C atoms
        atom_ids2 = []
        for r in res_list:
            atoms = r.fetch_atoms( atom_names_backup )
            for atom in atoms: atom_ids2.append( atom.id - 1 )
            
        if len(self.atom_ids)<4:
            random_int=rd.sample(range(0, len(atom_ids2)), 4-len(self.atom_ids))
            for i in random_int:
                self.atom_ids.append(atom_ids2[i])
    
        return self.atom_ids
    
    def fetch_residues(self):
        # get the residue_number of the ligand and the protein
        m=Model(self.input_pdb)
        for resi in m.residues:
            if resi.resname ==self.fragment_name:     
                self.lig_ind=resi.id
                self.n_res=self.lig_ind-1
    
        return self.n_res, self.lig_ind
    
    def input_contacts(self):
           
        n_res , lig_ind = self.fetch_residues()   
        ind_arr = np.zeros((n_res,2))
        for i in range(n_res):
            ind_arr[i][0] = lig_ind-1
            ind_arr[i][1] = i
        return ind_arr
    
    def get_coordinates(self):
        feat = coor.featurizer(self.input_pdb)
            
        if self.coordinates=='contacts':
            self.reduce_dimensions=True 
            feat.add_residue_mindist(residue_pairs=self.input_contacts(), scheme='closest-heavy', threshold=0.5)
        elif self.coordinates=='heteroatom':
            feat.add_selection(self.fetch_atoms())
        elif self.coordinates=='prot_frag':
            print 'to be implemented'
            # to be implemented
    
        self.inp = coor.source(self.traj_list, feat)
        X = self.inp.get_output()

        #reduce dimensions
        lag_time_tica=self.lagtime
        if self.reduce_dimensions:
            tica_obj = coor.tica(X, lag=lag_time_tica, kinetic_map=True)
            Y = tica_obj.get_output()  # get tica coordinates
        else:
            Y=X
        self.Y=Y    
        return Y
    
    def load_distances(self):
        feat2=coor.featurizer(self.input_pdb)
        feat2.add_residue_mindist(residue_pairs=self.input_contacts(), scheme='closest-heavy')
        P = coor.load(self.traj_list, feat2)
        
        return P
    
    def clustering(self):
        if self.cluster_method=='k_means':
            clustering = coor.cluster_kmeans(self.Y, k=self.n_clusters, max_iter=100, tolerance=1e-10)
            self.clustering_centers=clustering.clustercenters
            self.dtrajs=clustering.dtrajs
            n_clusters_local=self.n_clusters
        elif self.cluster_method=='uniform_time':
            clustering = coor.cluster_uniform_time(self.Y, k=self.n_clusters)
            self.clustering_centers=clustering.clustercenters
            self.dtrajs=clustering.dtrajs
            n_clusters_local=self.n_clusters
            
        if self.visual:
            self.plot_clustering(self.Y, clustering)
        
        return {'trajs':self.dtrajs, 'n_clusters':n_clusters_local, 'clustering_centers': self.clustering_centers, 'clustering': clustering}

    def plot_clustering(self, Y, clustering):
            ax, f = mplt.plot_free_energy(np.vstack(Y)[:,0], np.vstack(Y)[:,1])
            cc_x = clustering.clustercenters[:,0]
            cc_y = clustering.clustercenters[:,1]
            #ax.plot(cc_x,cc_y, linewidth=0, marker='o', markersize=5, color='black')
            #ax.set(xlabel='component 1', ylabel='component 2')
            return ax
            #ax.show()


class sort_cluster:
        
    '''
    In this class small clusters are removed and added to the next neighboring cluster. Then the distance between the 
    fragment and the protein is used to lump all unbound states together. 
    '''
    def __init__(self, clustering_class ):
        self.clustering_class = clustering_class
        self.lesscluster=0
        self.description = "In this class small clusters are removed and added to the next neighboring cluster. Then the distance between the fragment and the protein is used to lump all unbound states together. "
        
    def remove_small_clusters(self):
        small_clusters=[]
        for i in range(self.clustering_class.n_clusters):
            traj_lenght=0
            for n in range(len(self.clustering_class.Y)):
                #count how often the cluster appears in the list
                traj_lenght=traj_lenght+list(self.clustering_class.dtrajs[n]).count(i)
            if traj_lenght<5:
                small_clusters.append(i)
         
        all_states=[]
        for i in range(self.clustering_class.n_clusters+1):
            all_states.append(i)
            
        old_centers=self.clustering_class.clustering_centers
        cluster_list=np.setdiff1d(all_states,small_clusters)

        n_clusters_less=len(cluster_list)-1
        new_centers=[[] for m in range(n_clusters_less) ]
        n=0
        for i in range(self.clustering_class.n_clusters):
            if i in cluster_list:
                new_centers[n]=array(old_centers[i])
                n=n+1
        self.lesscluster=coor.assign_to_centers(data=self.clustering_class.Y, centers=new_centers)
        self.new_centers=new_centers
        self.n_clusters_less=n_clusters_less
    
    def assign_bound_unbound(self):
        Y=self.clustering_class.Y
        P=self.clustering_class.load_distances()
        i=0
        samples_water_2 = []
        watercluster_2=[[] for i in range(len(Y)) ]
        for n in range(len(Y)):
            for m in range(Y[n].shape[0]): 
                if min(P[n][m]) > 0.5:
                    # if distance of protein and fragment is greater the 0.5 nm, it is unbound
                    samples_water_2.append(len(self.new_centers))
                    m+=1
                else:
                    samples_water_2.append(self.lesscluster[n][m])
                    m+=1
            watercluster_2[n]=array(samples_water_2, dtype=int32)
            samples_water_2 = []
            n+=1
        return watercluster_2
    


class simple_MSM:
            
    '''
    In this class small clusters are removed and added to the next neighboring cluster. Then the distance between the 
    fragment and the protein is used to lump all unbound states together. 
    '''
    def __init__(self, clustering_class, cluster_assignment, lagtime, out_data, method='median', MSM_type='rev' ):
        self.clustering_class = clustering_class
        self.lesscluster=cluster_assignment
        self.distr=None
        self.lagtime=lagtime
        self.MSM_type=MSM_type
        self.method=method
        self.mm_distr=None
        self.out_data=out_data
        self.prob=0
        self.scoring=0
        self.score_states=[]
        self.simple_MSM_ranking={}
        self.description = "In this class small clusters are removed and added to the next neighboring cluster. Then the distance between the fragment and the protein is used to lump all unbound states together. "
        
    def create_MSM(self):
        if self.MSM_type=='rev':
            self.distr=msm.estimate_markov_model(self.lesscluster, lag=self.lagtime )
        elif self.MSM_type=='non_rev':
            self.distr=msm.estimate_markov_model(self.lesscluster, lag=self.lagtime, reversible=False )
        elif self.MSM_type=='bayes':
            self.distr=msm.bayesian_markov_model(self.lesscluster, lag=self.lagtime )
            
        self.mm_distr=self.distr.stationary_distribution
        
    def cluster_member_coor(self, timeseries):
        cluster_coor_list=[]
        for i in range(len(timeseries)):
            cluster_coor_list.append(self.clustering_class.Y[timeseries[i][0]][timeseries[i][1]])
        return cluster_coor_list
        
    def get_highest_populated(self):
        inp=self.clustering_class.inp
        self.structure_list=[]
        self.simple_MSM_ranking={}
        #get the sorting with the highest ranked first
        sorted_cluster=np.argsort(self.mm_distr)[::-1]
        #get the probability of the 5 highest clusters
        self.prob=np.sort(self.mm_distr)[::-1]
        self.scoring=self.prob
        
        # look for the best 5 cluster
        for n in range(5):
            key='rank_%i' % (n+1)
            rank_cluster=sorted_cluster[n]
            
            self.simple_MSM_ranking[key]=self.distr.active_state_indexes[rank_cluster]
            
           
            ### return as .xtc file 
            outfile = self.out_data+'/simple_msm_most_rank_'+str(n+1)+'.xtc'
            #coor.save_traj (inp,  cluster_list, outfile=outfile )
        
            #save all frames of small clusters
            if len(self.simple_MSM_ranking[key])<100:
                coor.save_traj(inp, self.simple_MSM_ranking[key], outfile=outfile )
                save_list=self.simple_MSM_ranking[key]
                #print 'short', save_list
                
            #save 100 random frames of large clusters
            else:
                save_list_rd=rd.sample(range(0, len(self.simple_MSM_ranking[key])), 100)
                save_list=[self.simple_MSM_ranking[key][i] for i in save_list_rd]
                #print save_list
                coor.save_traj (inp,  save_list, outfile=outfile ) 
                
        
            ###return representive structure as .pdb
            
            #get the according coordinates
            cluster_coor_list=self.cluster_member_coor(save_list)
            
            #calculate the mean positions
            if self.method=='average':
                mean_pos= np.average(cluster_coor_list, axis=0) # just average position
            elif self.method=='median':
                mean_pos = np.median(cluster_coor_list, axis=0) # just average position
            else:
                print 'Error: This method does not exists'
                
            #look for entry most closest to mean
            representative_index=min(range(len(cluster_coor_list)), key=lambda i: np.linalg.norm(cluster_coor_list[i]-mean_pos))
            representative=save_list[representative_index]
            #print cluster_coor_list[representative_index]
            outfile = self.out_data+'/simple_MSM_centroid_rank_'+str(n+1)+'.pdb'
            self.structure_list.append(outfile)
            self.score_states.append(outfile)
            coor.save_traj (inp, [representative], outfile=outfile )
            



class Scoring:
    '''
    A good binding pose is occupied for at least several ns. This class calculates how long the fragments stays
    in the bining pose an gives scores according to both the residence time and the MSM.
    '''
    
    def __init__(self, reference_quality, out_data, traj_list, gromacs_version=5 ):
        self.reference_quality = reference_quality
        self.traj_list=traj_list
        self.quality_list=[]
        self.gromacs_version=gromacs_version
        self.rank_quality_control=[]
        self.binding_dict={}
        self.out_data=out_data
       
     
    def visits_of_state(self, current_ref):
        
        max_binding_time=0
        binding_events=0
        non_bind_events=0
        unbind_events=0
        binding_lengh_list=[]
        
        #First calculate the rmsd to the cluster_centroid of each traj frame.
        rmsd_traj_list=Gromacs_commands(self.gromacs_version).calc_rmsd_traj(current_ref, self.traj_list, self.out_data)
        #Then look have often the binding site is found (rmsd < 0.3) and how long the fragment stays there
        for entry in rmsd_traj_list:
            
            f=open(entry,'r')
            content=f.readlines()
            rmsd_list=[]
            for n in range(len(content)):
                #the lines starting with # or @ have no rmsd information
                if content[n].startswith('#') or content[n].startswith('@'):
                    pass
                # split the other lines to get only the rmsd information and not the time
                else:
                    rmsd_list.append(float(str(content[n]).split()[1]))
            #now scan the list for binding events 
            if min(rmsd_list) > 0.3:
                pass
            
            else:
                i_start=0
                t=0
                while i_start < len(rmsd_list):
                    binding=False
                    i=i_start
                
                    while i <= len(rmsd_list)-5:
                        #for a binding event the fragment must at least stay 5 ns
                        if rmsd_list[i]<0.3 and rmsd_list[i+1]<0.3 and rmsd_list[i+2]<0.3 and rmsd_list[i+3]<0.3 and rmsd_list[i+4]<0.3:
                            binding_time=i
                            binding=True
                            binding_events+=1
                            m=binding_time
                            break
                        i=i+1
                        
                    if not binding:
                        #print "in trajectory %s there are no binding events." % (bound_traj)
                        non_bind_events+=1
                        break
                
                    while m >= binding_time and m < len(rmsd_list):
                        if rmsd_list[m]> 0.5:
                            last_bound=m
                            unbind_events+=1
                            break
                        m=m+1
                        if m==len(rmsd_list):
                            last_bound=m
                            'print reached last'
    
                    binding_lenght=last_bound-binding_time
                    binding_lengh_list.append(binding_lenght)
                    
                    #look when are start and end binding time for longest bound state
                    #for the B-factor analysis
                    if binding_lenght> max_binding_time:
                        binding_lenght=max_binding_time
                        start=binding_time
                        end=last_bound
                        traj=self.traj_list[rmsd_traj_list.index(entry)]
            
                    i_start=last_bound
                
            #remove the old rmsd files
            subs = sbp.Popen(['rm', entry])
            #write start and end time of longest traj to a dict
        try:
            self.binding_dict[current_ref]=[traj, start, end]
        except:
            self.binding_dict[current_ref]=[None, 0, 0]
        
        return binding_events, unbind_events, non_bind_events, binding_lengh_list
        
    def calc_score(self, current_ref):
        binding_events, unbind_events, non_bind_events, binding_lengh_list=self.visits_of_state(current_ref)
        try:
            accum_binding_time= np.cumsum(binding_lengh_list)[-1]
        except:
            accum_binding_time=0
        #binding traj gives a bonus of 1, unbinding a bonus of 0.1. Traj that do not bind give a malus of -0.5. 
        #This is then weighted by the total binding time
        rank_score=(1*binding_events+0.1*unbind_events-0.5*non_bind_events)*accum_binding_time      
        if rank_score<0:
            rank_score=0
        return rank_score
            
    def quality_score(self):
        for entry in self.reference_quality:
            self.quality_list.append(self.calc_score(entry))
            
        if np.cumsum(self.quality_list)[-1]==0:   # do not devide by 0
            self.rank_quality_control=self.quality_list
        else: #weighted scoring
            self.rank_quality_control=(self.quality_list/np.cumsum(self.quality_list)[-1])*100
            
    def plot_score(self):
        print 'test'
            


class Kinetic_clustering:             
    '''
    This class aimes to lump similar bound states together to get a better model. 
    '''
    def __init__(self, Markov_class, out_data, nstates_cg=10, n_corestates=5, gromacs_version=5, max_trace_runs=250000, PCCA_exclusion_value=0.0, max_lenght_factor=0.4 ):
        self.m=Markov_class
        self.max_lenght_factor=max_lenght_factor
        self.remove_unlikely=True
        self.PCCA_exclusion_value=PCCA_exclusion_value
        self.max_trace_runs=max_trace_runs
        self.nstates_cg=nstates_cg
        self.lumpstates=[]
        self.method='median'
        self.structure_list=[]
        self.n_corestates=n_corestates
        self.rmsd_cg=0
        self.out_data=out_data
        self.gromacs_version=gromacs_version
        self.description = "In this class small clusters are removed and added to the next neighboring cluster. Then the distance between the fragment and the protein is used to lump all unbound states together. "
        
    def PCCA_states(self):
        self.m.distr.pcca(self.nstates_cg)
        pcca_dist=self.m.distr.metastable_distributions
        membership=self.m.distr.metastable_memberships
        sets=self.m.distr.metastable_sets
        
        # exclude clusters that have a low probability to be in the macrostate
        if self.remove_unlikely==True:
            sets2=[]
            for l in range(len(sets)):
                test_set=[]
                for m in range(len(sets[l])):
                    if membership[m][l]> self.PCCA_exclusion_value:
                        test_set.append(sets[l][m])
                sets2.append(test_set)
                
            self.kinetic_states=sets2 
        else:
            self.kinetic_states=sets 

        
    def maxtrace_states(self):
        
        ''' 
        code from Roberto Covino
        Input is the microscopic propagator M, number of reduced states N,
         NC number of core states and number of Monte Carlo sweeps

         Add arbitrary number of CG states.
         Hande expection of singular Y
         '''
        from scipy.stats import uniform
        M=self.m.distr.transition_matrix
        N=self.nstates_cg 
        NC=self.n_corestates 
        nmc=self.max_trace_runs
        
        temp0 = np.float( 10**-2 )          # starting temperature
        temp1 = np.float( 10**-10 )         # target (final) temperature
        alpha = np.log( temp1 / temp0 ) / nmc
        step_array = np.arange( 1, nmc + 1 )
        temp = temp0 * np.exp( alpha * step_array ) # decrease of effective temperature


        n = M.shape[ 0 ]

        ''' Silly Billy optimization '''
        idx = np.random.randint( N, size=n )
        A = np.identity( N, dtype=int )[ idx ] # randomly initialize A
        Ainit = np.copy( A )

        maxTrax = np.zeros( nmc + 1 )

        T = np.dot( np.transpose(A), np.dot( M, A) )   # reduced transition matrix
        T = T / T.sum( axis=0, dtype=float )
        # calculate and sort eigenvalues

        maxTrax[ 0 ] = np.sum( np.sort(np.real(eigvals( T )))[::-1][ :NC ] )

        '''
        vals, vecs = eig( M )
        peq = np.real( vecs[:,0] ) / np.sum( np.real( vecs[:,0] ) )
        # the eigenvectors returned by numpy are not normalized! You should 
        sort the eivals and eigec the same way!
        Dn = peq.reshape((-1,1)) * np.identity( n )   # Matrix with peq 
        elements on the diag
        X = np.identity( n ) - M + np.ones(( n, n )) * peq.reshape((-1,1))
        Xinv = inv( X )

        while True:  # tries different initializations until matrix is not 
        singular

        idx = np.random.randint( N, size=n )
        A = np.identity( N, dtype=int )[ idx ] # randomly initialize A

        V = np.dot( Xinv, Dn )
        Y = np.dot( np.transpose(A), np.dot( V, A) )
        try:
            Yinv = inv(Y)
        except:
            continue
        break

        Ainit = np.copy( A )

        Peq = np.dot( peq, A )
        DN = Peq.reshape((-1,1)) * np.identity( N )

        Z = np.dot( DN, Yinv )

        maxTrax = np.zeros( nmc + 1 )
        maxTrax[ 0 ] = N + 1 - np.trace( Z )

        maxT = np.identity( N ) + np.ones(( N, N )) * Peq.reshape((-1,1)) - Z
     '''
        ru = uniform.rvs( size = nmc )
        w = np.zeros( nmc )

        for step in step_array:


            idx_flip = randint( n, size=1 )[0]
            trialA = np.copy( A )

            trialA[ idx_flip ] = permutation( trialA[ idx_flip ] ) # add a check to try untile
                                                                 # the permutation is different

            T = np.dot( np.transpose(trialA), np.dot( M, trialA) )
            T = T / T.sum( axis=0, dtype=float )
            vals = np.sort(np.real(eigvals( T )))[::-1]
            #print vals[ :NC ]
            trialTrax = np.sum( vals[ :NC ] )

            ''' Simulated annealing '''
            w[ step - 1 ] = np.exp(( trialTrax - maxTrax[ step - 1 ] ) / temp[ step - 1 ] )

            if w[ step - 1 ] > ru[ step - 1 ]:                        # accept the move
                maxTrax[ step ] = trialTrax
                A = np.copy( trialA )
            else:
                maxTrax[ step ] = maxTrax[ step - 1 ]
                
        cluster_dict=[]
        for i in range(len(A[0])):
            cluster_list=np.where( A[:,i] == 1)[0]
            cluster_dict.append(cluster_list)
        self.kinetic_states=cluster_dict

        
    def lump_cluster_all(self):
        #Chose the macrocluster to lump. Very large clusters will not be lumped.
        self.lumpstates=[]
        for k in range(len(self.kinetic_states)):
            if (1<len(self.kinetic_states[k])<=self.max_lenght_factor*self.m.clustering_class.n_clusters) :
                self.lumpstates.append(self.kinetic_states[k])
        
    def lump_cluster_max(self):
        #Chose the macrocluster to lump. Only the macrocluster with high ranked structures inside 
        # are considerd. Very large clusters will not be lumped
        self.lumpstates=[]
        sets_highest_prob=[]
        highest_cluster=np.argsort(self.m.mm_distr)[::-1][:10]
        
        #check which macrocluster contain the highest ranked ones 
        for entry in highest_cluster:
            for n in range(len(self.kinetic_states)):
                if entry in self.kinetic_states[n]:
                    sets_highest_prob.append(n)
        sets_analysis=list(set(sets_highest_prob))
                        
        for k in range(len(self.kinetic_states)):
            if k in sets_analysis:
                self.lumpstates.append(self.kinetic_states[k])
        
    def lump_cluster_max_small(self):
        #Chose the macrocluster to lump. Only the macrocluster with high ranked structures inside 
        #and small cluster_size are considerd.
        self.lumpstates=[]
        sets_highest_prob=[]
        highest_cluster=np.argsort(self.m.mm_distr)[::-1][:10]
        
        #check which macrocluster contain the highest ranked ones 
        for entry in highest_cluster:
            for n in range(len(self.kinetic_states)):
                if entry in self.kinetic_states[n]:
                    sets_highest_prob.append(n)
        sets_analysis=list(set(sets_highest_prob))
                        
        for k in range(len(self.kinetic_states)):
            if k in sets_analysis and (1<len(self.kinetic_states[k])<=self.max_lenght_factor*self.m.clustering_class.n_clusters):
                self.lumpstates.append(self.kinetic_states[k])
                
    def cluster_member_coor(self, timeseries):
        cluster_coor_list=[]
        for i in range(len(timeseries)):
            cluster_coor_list.append(self.m.clustering_class.Y[timeseries[i][0]][timeseries[i][1]])
        return cluster_coor_list
    
    
    def k_means_cg_nodes(self, reference_list, rmsd_cg):
        # implement something like k_mean clustering 
        #select a random startin node
        clustercenter=[]
        start_node=np.random.randint(0, len(reference_list) )
        clustercenter.append(start_node)

        #identify the first cluster farest away
        highest_center_old= max(rmsd_cg[start_node])
        cluster_assignment_list=[start_node for n in range(len(reference_list))]
        if highest_center_old >0.3:
            new_node=argmax(rmsd_cg[start_node])

        while highest_center_old > 0.3 :     
    
            clustercenter.append(new_node)
    
            # reassign the clusters 
            cluster_assignment_list=[]
            for i in range(len(reference_list)):
                distance_new=[]
                for j in clustercenter:
                    distance_new.append(rmsd_cg[i][j])
                assignment=argmin(distance_new)
                cluster_assignment=clustercenter[assignment]
                cluster_assignment_list.append(cluster_assignment)
    
    
            #identify again largest cluster 
            #just compare clusters within the same macrocluster:
        
            highest_center_old=0

            for i in clustercenter:
                highest_distance=[]
                for j in range(len(reference_list)):
                    if cluster_assignment_list[j]==i :
                        highest_distance.append(rmsd_cg[i][j])
                highest_center=max(highest_distance)


                if highest_center>highest_center_old:
                    highest_center_old=highest_center
                    new_node_pos = [b for b, x in enumerate(rmsd_cg[i]) if x == highest_center_old]
                    free_nodes=list(set(new_node_pos)-set(clustercenter))
                    #first node will always be chacked by default, however after subtracting the clustercenter
                    #the list might be empty and cause an error, so introduce a pass!
                    try:
                        new_node=free_nodes[0]
                    except:
                        pass   

        return cluster_assignment_list 

    
        
    def similarity_restraint(self, method='median'):
        #check the structural similarity between the clusters. Remove structural outliers.
        self.method=method
        similarity_list=[]
        inp=self.m.clustering_class.inp
        lumping_list=self.lumpstates
        #do the similarity restrain for all the sublists
        for sub_lumping_list in lumping_list:
            self.structure_list=[]
            for cluster in sub_lumping_list:
                
            #create centroid of each cluster
            
                lumpin_coor=self.m.distr.active_state_indexes[cluster]
            
                ### return as .xtc file 
                outfile = self.out_data+'/kinetic_clustering_'+str(cluster)+'.xtc'
        
                #save all frames of small clusters
                if len(lumpin_coor)<100:
                    #coor.save_traj(inp, lumpin_coor, outfile=outfile )
                    save_list=lumpin_coor
                #print 'short', save_list
                
                #save 100 random frames of large clusters
                else:
                    save_list_rd=rd.sample(range(0, len(lumpin_coor)), 100)
                    save_list=[lumpin_coor[i] for i in save_list_rd]
                    #coor.save_traj (inp,  save_list, outfile=outfile ) 
            
                #return representive structure as .pdb
            
                #get the according coordinates
                cluster_coor_list=self.cluster_member_coor(save_list)
            
                #calculate the mean positions
                if self.method=='average':
                    mean_pos= np.average(cluster_coor_list, axis=0) # just average position
                elif self.method=='median':
                    mean_pos = np.median(cluster_coor_list, axis=0) # just average position
                else:
                    print 'Error: This method does not exists'
                
                #look for entry most closest to mean
                representative_index=min(range(len(cluster_coor_list)), key=lambda i: np.linalg.norm(cluster_coor_list[i]-mean_pos))
            
                representative=save_list[representative_index]
                #print cluster_coor_list[representative_index]
                outfile = self.out_data+'/kinetic_clustering_centroid_'+str(cluster)+'.pdb'
                self.structure_list.append(outfile)
                coor.save_traj (inp, [representative], outfile=outfile )
            
            #calculate the rmsd between the structures
            rmsd_cg=np.zeros((len(sub_lumping_list), len(sub_lumping_list)))

            # use only the half matrix cause 1.2=2.1
            for n in range(len(sub_lumping_list)):
                m=0
                while m <= n:
                    rmsd_cg[n][m]=Gromacs_commands(self.gromacs_version).calc_rmsd_ref(self.structure_list[n], self.structure_list[m], self.out_data)
                    rmsd_cg[m][n]=rmsd_cg[n][m]
                    #remove the old rmsd files
                    m+=1
            self.rmsd_cg=rmsd_cg  
            #exclude clusters that are not structural similar enough
    
            cluster_assignment_list=self.k_means_cg_nodes(self.structure_list, rmsd_cg)
            
            #new list for lumping states together 
            for m in set(cluster_assignment_list):
                small_list=[]
                indices = [i for i, x in enumerate(cluster_assignment_list) if x == m]
                for k in indices:
                    small_list.append(sub_lumping_list[k])
                #single clusters are not lumped
                if len(small_list)>1:
                    similarity_list.append(small_list) 
            
        self.lumpstates= similarity_list

        
    def lump_states(self):
        #Create new assignments with the lumped states
        # do this for all states that should be lumped together
        
        #msm buidling does not work with just the active trajectories, 
        # so use the full one 
        traj_copy=self.m.distr.discrete_trajectories_full
        #update the lumpstates to full traj
        update_lumpstates=[ self.m.distr.active_set[x] for x in self.lumpstates]
    
        for l in range(len(update_lumpstates)):
            current_lump=update_lumpstates[l]
    
            #uptdate lesscluster if possible
            if l>0:
                traj_copy=trajecory_cg
        
            i=0
            cg_list = []
            cg_traj=[[] for i in range(len(traj_copy)) ]
            for n in range(len(traj_copy)):
                for m in range(len(traj_copy[n])): 
        
                    if traj_copy[n][m] in  current_lump:
                    # g is necessary if lumping is done in multiple macrostates 
                        #new, lumped clusters get new, (larger), clusternumbers
                        b = self.m.distr.nstates_full+1+l
                        cg_list.append(b)
                    else:
                        c=traj_copy[n][m]
                        cg_list.append(c)
                cg_traj[n]=array(cg_list, dtype=int32)
                cg_list = []
            trajecory_cg=cg_traj
        
        return trajecory_cg



class output_section:
    '''
    This class can give out the plot of the score and the B-factor of the ligand.
    '''
    def __init__(self, out_data, input_pdb, gromacs_version=5):
        self.B=False
        self.out_data=out_data
        self.scoring=0
        self.input_pdb=input_pdb
        self.gromacs_version=gromacs_version
        self.score_states=[]
        self.description= 'This class can give out the plot of the score and the B-factor of the ligand.'
        
    def B_factor(self, quality_score_class):
        # create B-factor for bound states 
        self.B=True
        b_fac_dict=quality_score_class.binding_dict
        for rank in range(1,6):
            name=  self.out_data+'/simple_MSM_centroid_rank_%i.pdb' % (rank)     
            b_list=b_fac_dict[name]
            #times have to be multiplied by 1000 to get the ps unit
            traj, begin_time, end_time= b_list[0], b_list[1]*1000, b_list[2]*1000
            if traj!= None:
                Gromacs_commands(self.gromacs_version).calc_b_factor(traj, self.input_pdb, rank, begin_time, end_time, self.out_data)
            else:
                print 'No binding trajectories detected for cluster rank %i' % (rank)
        
    def plot_ranking(self, quality_rank, MSM_class):
        score_states_prev=[]
        #calculate MSM_ranking first
        MSM_score=100*MSM_class.prob[:5]/np.cumsum(MSM_class.prob[:5])[-1]
        
        #add to quality score
        scoring = quality_rank+MSM_score
        sort_scoring=np.sort(scoring)[::-1]
        self.scoring=sort_scoring
        
        #plot score
        pylab.plot( sort_scoring, range(1,6), 'o' )
        pylab.legend(bbox_to_anchor=(2.,1))
        pylab.grid(lw=2)
        pylab.xlabel('score', fontsize=16)
        pylab.ylabel('rank', fontsize=16)
        pylab.ylim(0,5)
        pylab.show()
        
        #print representative 
        n=0
        for item in MSM_class.structure_list:
            sort_scoring_list=list(sort_scoring)
            scoring_rank=sort_scoring_list.index(scoring[n])
            copy_name='%s/scoring_rank%i.pdb' % ( self.out_data, scoring_rank+1)
            subs = sbp.Popen(['cp', item, copy_name])
            score_states_prev.append(copy_name)
            if self.B:
                #also copy the B_factor
                item_B='%s/B_fac_%i.pdb' % ( self.out_data, n+1)
                copy_name_B='%s/B_fac_scoring_rank%i.pdb' % ( self.out_data, scoring_rank+1)
                subs = sbp.Popen(['cp', item_B, copy_name_B])
            n+=1
        self.score_states=list(sort(score_states_prev))



class error_estimation:
    '''

    The error sensitivity analysis can be used for adaptive sampling.

    '''

    def __init__(self, out_data, method='median'):
        self.sensitivity_matrix=0
        self.state_errors=0
        self.out_data=out_data
        self.error_cluster=0
        self.Bays_MSM_class=None
        self.structure_list=[]
        self.error_ranking={}
        self.method=method
        self.description = "The error sensitivity analysis can be used for adaptive sampling."
        
    def sensitivity_analysis(self, Bays_MSM_class,h=0):
        # here we do the sensitivity_analysis for the h-th eigenvector (default h=0, equilibrium) 
        # senitivity of the eigenvectector will be included later
        
        #the sensitivity analysis will us tell how to weight the clusters
        self.Bays_MSM_class=Bays_MSM_class
        P=Bays_MSM_class.distr.transition_matrix
        self.sensitivity_matrix=self.get_eigenvalue_sensitivity(P,h)
        
        #an Bayesian MSM will be used to define the errors 
        self.state_errors=Bays_MSM_class.distr.sample_std('stationary_distribution')

        # we want to resample highly important states with high errors, so we
        # multiply the errors with the weight and resample the highest states 
        
        #sum over the columns to get the sensitivity of one cluster
        k=len(P)
        self.sensitivity_cluster=np.zeros(k)
        for i in range(k):
            sensi=0
            for j in range(k):
                    sensi=sensi+self.sensitivity_matrix[i][j]
            self.sensitivity_cluster[i]=sensi
    
        self.error_cluster=np.zeros(k)
        for item in range(len(self.sensitivity_cluster)):
            self.error_cluster[item]=self.sensitivity_cluster[item]*self.state_errors[item]
            
        # sort the clusters after their error and write out the highest 5 error states 
        self.sort_error=np.sort(self.error_cluster)[::-1]
        self.sort_error_states=self.error_cluster.argsort()[::-1]
        # write out the states 
        self.get_highest_error(self.sort_error_states[:5])
            
 
    def get_highest_error(self, cluster_list):  
        #write out states with highest error
        inp=self.Bays_MSM_class.clustering_class.inp

        for n in range(5):
            key='error_%i' % (n+1)
            rank_cluster=cluster_list[n]
            
            self.error_ranking[key]=self.Bays_MSM_class.distr.active_state_indexes[rank_cluster]
            
           
            ### return as .xtc file 
            outfile = self.out_data+'/error_rank_'+str(n+1)+'.xtc'
        
            #save all frames of small clusters
            if len(self.error_ranking[key])<100:
                coor.save_traj(inp, self.error_ranking[key], outfile=outfile )
                save_list=self.error_ranking[key]
                
            #save 100 random frames of large clusters
            else:
                save_list_rd=rd.sample(range(0, len(self.error_ranking[key])), 100)
                save_list=[self.error_ranking[key][i] for i in save_list_rd]
                coor.save_traj (inp,  save_list, outfile=outfile ) 
                
            ###return representive structure as .pdb
            
            #get the according coordinates
            cluster_coor_list=self.cluster_member_coor(save_list)
            
            #calculate the mean positions
            if self.method=='average':
                mean_pos= np.average(cluster_coor_list, axis=0) # just average position
            elif self.method=='median':
                mean_pos = np.median(cluster_coor_list, axis=0) # just average position
            else:
                print 'Error: This method does not exists'
                
            #look for entry most closest to mean
            representative_index=min(range(len(cluster_coor_list)), key=lambda i: np.linalg.norm(cluster_coor_list[i]-mean_pos))
            representative=save_list[representative_index]
            outfile = self.out_data+'/error_centroid_rank_'+str(n+1)+'.pdb'
            self.structure_list.append(outfile)
            coor.save_traj (inp, [representative], outfile=outfile )  
            
    def cluster_member_coor(self, timeseries):
        cluster_coor_list=[]
        for i in range(len(timeseries)):
            cluster_coor_list.append(self.Bays_MSM_class.clustering_class.Y[timeseries[i][0]][timeseries[i][1]])
        return cluster_coor_list
        
    def get_eigenvalue_sensitivity(self,P,h):
        # P is the probability matrix, you get the sensitivity matrix for the h-th eigenvalue 
        
        #calculate Eigenvalues and Eigenvectors 
        eigenval, eigenvec=np.linalg.eig(P)
        # define matrix A with A=P-eigenval*I
        k=len(P)
        A=P-eigenval[h]*np.identity(k)
    
        #check the rank of the Matrix
        if np.linalg.matrix_rank(A)!= k-1:
            print 'error: rank of the matrix is not the lenght of the matrix minus 1'
    
        # for the sensitivity analysis according to N.S. Hinrichs and V.S Pande 2007 we want to have the unity values
        #along the upper and the zero elemts along the lower triangular matrix
        #Therefore we first have to calculate L and U for the transposed matrix B and then transpose them
        B=A.transpose()

        mut, L_B, U_B=scipy.linalg.lu(B)

        # now transpose them
        L_A=U_B.transpose()
        U_A= L_B.transpose()
    
        # now we have to solve the following equations:
        # U_A*x=ek where ek is the column vector corresponding to the kth column of the identity matrix
        # L_A.transpose()*xa=0 with setting the kth element of xa to 1 to get a nontrivial result 

        #equation 1
        ek=np.reshape(np.identity(k)[k-1], (k, 1))
        x=np.linalg.solve(U_A, ek)

        #equation 2

        #l_kk is zero within the numerical errors
        L_A[k-1,k-1] = 0

        # the last row of L_A.transpose() contains only zeros, so xa_k can be chosen randomly. We define it as 1.
        # To solve the equation we move the values containing xa_k to the righthand site and the solve the reduced matrix

        #solve the equation L_A.transpose()*xa=b with reduced matrix
        #define b
        b=np.zeros(k-1)
        for n in range(k-1):
            b[n]=-L_A.transpose()[n][k-1]

        #reduce the matrix 
        L_A_reduced= L_A.transpose()[np.ix_(range(k-1),range(k-1))]

        #solve
        xa_prev=np.linalg.solve(L_A_reduced, b)
        xa=np.append(xa_prev, [1]).real  #change
    
        # calculate the eigenvalue sensitivity with equation A12 (paper)
        normalization=np.dot(xa.transpose(), x).real
        sensitivity_matrix=np.zeros([k,k])
        for i in range(k):
            for j in range(k):
                sensitivity_matrix[i][j]=xa[i]*x[j].real/normalization
            
        return sensitivity_matrix

    def get_eigenvector_sensitivity(P,h,i,j):
        # P is the probability matrix, you get the sensitivity for the h-th eigenvector
        #on the transition for i to j
        
        sensitivity_matrix=self.get_eigenvalue_sensitivity(P,h)
        #calculate Eigenvalues and Eigenvectors 
        eigenval, eigenvec=np.linalg.eig(P)
        # define matrix A with A=P-eigenval*I
        k=len(P)
        A=P-eigenval[h]*np.identity(k)
    
        #update A to solve the equations 
        newrow =eigenvec[h]
        A_new= np.vstack([A, newrow])
    
        #equation1
        d=np.linalg.lstsq(A_new, eigen_new)
        second=(1/sensitivity_matrix[i][j])*d[0]
    
        #equation2
    
        ei=np.reshape(np.identity(k)[i-1], (k, 1))
        ei_new=np.vstack([ei, 0])
        ci=np.linalg.lstsq(A_new, -ei_new)
        first=(ci[0]/eigenvec[0,j])
    
        eigenvector_sensitivity=first+second
    
        return eigenvector_sensitivity
    


class integration_and_bootstrapping:
    '''
    In this class the nuisance parameters like lagtime or the number of clusters is integrated out.  Also it is 
    possible to estime the error due to the limited simulation time by bootstrapping. 
    '''

    def __init__(self, out_data, gromacs_version=5, reference=None):
        self.description = "In this class the nuisence parameters like lagtime or the number of clusters is integrated out.  Also it is possible to estime the error due to the limited simulation time by bootstrapping."
        self.bootstrap_dict={}
        self.bootstrap_cluster_dict={}
        self.bootstrap_rank_dict={}
        self.bootstrap_rmsd_dict={}
        self.runnumber=0
        self.reference=reference
        self.out_data=out_data
        self.MSM_ranking=None
        self.gromacs_version=gromacs_version
        self.reference_folder=out_data+'/bootstrap_references'
        
        
    def compare_top_results(self, MSM_ranking, runnumber):
	self.runnumber=runnumber
        self.MSM_ranking=MSM_ranking
        # make dictionary with representative clusters and accumulate their score/prob
        i=0
        score_states=self.MSM_ranking.score_states
        scoring_perc=self.MSM_ranking.scoring
        for name in score_states:
            #calculate rmsd between the existing 
            ref1=name
            if len(self.bootstrap_dict)>0:
                included=False
                for key in self.bootstrap_cluster_dict:                                

                    ref2=key
                    
                    boot_rmsd=Gromacs_commands(self.gromacs_version).calc_rmsd_ref(ref1, ref2, self.out_data)
                    if boot_rmsd<0.3:
                          
                        self.bootstrap_dict[key].append(scoring_perc[i])
                        self.bootstrap_rank_dict[key].append(i+1)                       
                        self.bootstrap_cluster_dict[key].append(self.runnumber)
                        included=True
                        try:
                            self.bootstrap_rmsd_dict[key].append(Gromacs_commands(self.gromacs_version).calc_rmsd_ref(ref1, self.reference, self.out_data))
                        except:
                            pass
                
                                           
                if included==False:
                    key_new='%s/bootstrap%i.pdb' % (self.reference_folder, len(self.bootstrap_dict)+1)

                    self.bootstrap_dict[key_new]=[scoring_perc[i]]
                    self.bootstrap_cluster_dict[key_new]=[self.runnumber] 
                    self.bootstrap_rank_dict[key_new]=[i+1]
                    try:
                        self.bootstrap_rmsd_dict[key_new]=[Gromacs_commands(self.gromacs_version).calc_rmsd_ref(ref1, self.reference, self.out_data)]
                    except:
                        pass
                                            
                    # copy reference into the folder 
                    command_copy='cp %s %s' % (name, key_new )
                    os.system(command_copy)
                                                

            else:
                key_1='%s/bootstrap1.pdb' % (self.reference_folder)
                
                self.bootstrap_dict[key_1]=[scoring_perc[i]]
                self.bootstrap_cluster_dict[key_1]=[self.runnumber]
                self.bootstrap_rank_dict[key_1]=[i+1]
                #do also an rmds analysis, if a reference is defined 
                try:
                    self.bootstrap_rmsd_dict[key_1]=[Gromacs_commands(self.gromacs_version).calc_rmsd_ref(ref1, self.reference, self.out_data)]
                except:
                    pass
                
                #create_directory
                if not os.path.exists(self.reference_folder):
                    os.makedirs(self.reference_folder)
                # copy reference into the folder 
                command_copy='cp %s %s' % (name, key_1)
                os.system(command_copy)
            i=i+1

    def integrate_nuisance(self, n_cluster_pos,lagtime_pos, traj_list):
        print 'Test'
    
    def bootstrap_simulation_data(self):
        print 'Test'

    def bootstrap_mean_std(self, runnumber):

	self.row_mean=[]
	self.row_std=[]
	self.row_std_err=[]
	self.row_rmsd_mean=[]
	self.row_rmsd_std=[]
	self.row_rmsd_std_err=[]

	runs=runnumber-1


	for key in self.bootstrap_dict:
    		row=[]
    		row2=[]
    		for k in range(1, runs+1):
        		# the variable var accumulates all scores that belong to the same state (e.g state was split up)
        		var=0
        		for l in range(len(self.bootstrap_dict[key])):
            			if self.bootstrap_cluster_dict[key][l] ==k:
                			var=var+self.bootstrap_dict[key][l]
        		row.append(var)

    		#save the mean score
    		self.row_mean.append(np.mean(row))
    		#save the std
    		self.row_std.append(np.std(row))
    		#save the std-err
    		self.row_std_err.append(np.std(row)/np.sqrt(runs))

		if not (self.reference is None): 

    			#save the mean rmsd 
    			self.row_rmsd_mean.append(np.mean(self.bootstrap_rmsd_dict[key])) 
    			#calculate std of rmsd 
    			self.row_rmsd_std.append(np.std(self.bootstrap_rmsd_dict[key]))
    			#calculate stderr of rmsd 
    			self.row_rmsd_std_err.append(np.std(self.bootstrap_rmsd_dict[key])/np.sqrt(runs))
		else:
    			#save the mean rmsd 
    			self.row_rmsd_mean.append(1) 
    			#calculate std of rmsd 
    			self.row_rmsd_std.append(0)
    			#calculate stderr of rmsd 
    			self.row_rmsd_std_err.append(0)

    def bootstrap_plot_results(self, modus='std_err'):

	if modus == 'std_err':
    		xerr=self.row_std_err
    		yerr=self.row_rmsd_std_err
	elif modus =='std':
    		xerr=self.row_std
    		yerr=self.row_rmsd_std
    
	pylab.errorbar(self.row_mean, self.row_rmsd_mean, xerr=xerr, yerr=yerr, fmt='o')
	pylab.xlabel('score')
	pylab.ylabel('mean rmsd [nm]')

        





