import numpy as np
from read_kevin_bv_tracks import *
import pdb
"""
Function for reading clusters already produced and finding the BVs that were included in the clustering
Function for plotting cluster densities
Function for doing Kmeans clustering
"""

# some constants related to clustering
K=3
nsamples=6  # number of location samples from each BV used for clustering
lon_ix=np.arange(nsamples)*2
lat_ix=np.arange(nsamples)*2+1
# when we did the clustering we looked at locations in the region defined by lon_range and lat_range
# a few BVs go off in strange directions towards the end of their life so we don't include that part
cluster_lon_range=[90,120]
cluster_lat_range=[-10,20]
cluster_dir='/gws/nopw/j04/forsea/users/jcrook/BorneoVortex/Clusters/' # where we store the output

"""
read_clusters()
This reads the data we produced for kmeans clustering for k=3 plus sub k=3 and given a list of BVs checks which of these were involved in the clustering
This data consists of:
cluster_labels - an array which is the same size as the True BVs (i.e. init_in_reg and with no TCs)
cluster_names - an array holding the names of the clusters so they are more meaningful (eg. S
cluster_centers - an array holding the 6 centroids for each cluster [nbvs, nsamples*2] where lon_ix are the indices
                 in which the longitudes are stored and lat_ix are the indices where the latitudes are stored
cluster_bvids - ids of the BVs that were used to do the clustering
cluster_bv_startdates - the startdates of the BVs that were used to do the clustering

These latter two allow us to check that the cluster_labels match the BVs

Inputs:
    bv_tracks850 is the database of BVs from which this finds those that were clustered
    cluster_dir is where to read the data from and defaults to the cluster_dir defined above
Returns:

"""
def read_clusters(bv_tracks850, cluster_dir=cluster_dir):

    print('reading cluster data from', cluster_dir)
    cluster_labels=np.load(cluster_dir+'kmeans3_plus_sub.npy',allow_pickle=True) # the assigmnent of cluster to the BVs
    cluster_bvids=np.load(cluster_dir+'kmeans_bvids.npy',allow_pickle=True) # the BV track_ids
    cluster_bv_startdates=np.load(cluster_dir+'kmeans_bv_startdates.npy',allow_pickle=True) # the BV startdates
    cluster_centers=np.load(cluster_dir+'kmeans3_plus_sub_centers.npy', allow_pickle=True) # the cluster centroids

    # use the centers to determine the names
    max_cluster=np.amax(cluster_labels)
    cluster_names=['']*K*2
    letters=['a','b','c']
    # get the cluster which was sub-clustered
    unique=np.unique(cluster_labels)
    print(unique)
    sub_cix=np.where(np.asarray([cix not in unique for cix in np.arange(K*2)]))[0][0]
    print('cluster that was subclustered is', sub_cix)

    for c in range(max_cluster+1):
        lons=cluster_centers[c,lon_ix]
        lats=cluster_centers[c,lat_ix]
        if c==sub_cix:
            # this is the cluster that was subclustered
            cluster_names[c]='Cluster {c:d}'.format(c=c+1)
        elif np.amin(lats)>5:
            if lons[-1]<100:
                cluster_names[c]='S China Sea'
            else:
                cluster_names[c]='N Borneo'
        else:
            if np.amin(lats)<0:
                cluster_names[c]='C-shaped'
            else:
               if np.amax(lons)>110:
                    cluster_names[c]='NW Borneo'
               else:
                    cluster_names[c]='W Borneo'

    # now work out which of the given BVs were clustered
    nbvs=len(bv_tracks850)
    bv_clustered=np.zeros(nbvs,bool)
    next_cl_bv=0
    for b in range(nbvs):
        if bv_tracks850[b].track_id==cluster_bvids[next_cl_bv] and bv_tracks850[b].track_times[0]==cluster_bv_startdates[next_cl_bv]: 
            bv_clustered[b]=True
            next_cl_bv+=1
        else:
            bv_clustered[b]=False
    ix=np.where(bv_clustered)
    clust_bv_tracks=bv_tracks850[ix]
    if len(clust_bv_tracks)!=len(cluster_labels):
        print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
        print("WARNING: number of clustered bvs does not match number of cluster labels!!!!!!!!!!")
        
    return cluster_labels, cluster_centers, np.asarray(cluster_names), sub_cix, clust_bv_tracks
"""
plot_cluster_densities()
function to plot the track densities of the clusters - creates the figure for the user to save or display
Inputs:
    bv_tracks - the list of BVs that were clustered
    cluster_labels, cluster_centers, cluster_names, sub_cix - as returned from read_clusters
"""
from get_track_density import *
def plot_cluster_densities(bv_tracks, cluster_labels, cluster_centers, cluster_names, sub_cix):
    
    nbvs=len(bv_tracks)
    # work out the number of days over which the BVs exist - needed for densities
    start_dates=np.asarray([bv.track_times[0] for bv in bv_tracks])
    end_dates=np.asarray([bv.track_times[-1] for bv in bv_tracks])
    start_date=np.amin(start_dates)
    end_date=np.amax(end_dates)
    deltat=(end_date-start_date)
    ndays=deltat.days # this is including all days of year but we only look at BVs in Oct-Mar
    print(ndays, start_date.strftime('days from %d/%m/%Y '), end_date.strftime('to %d/%m/%Y'))
    # work out how many just in Oct-Mar season
    all_days=np.asarray([start_date+dt.timedelta(days=day) for day in range(ndays)])
    in_seas=np.asarray([(date.month<4) | (date.month>=10) for date in all_days])
    ix_seas=np.where(in_seas)
    ndays=len(ix_seas[0])
    print(ndays, 'days in Oct-Mar')

    # create grid on which to get clusters
    nlon=(cluster_lon_range[1]+5-cluster_lon_range[0])+1 # get densities to 125E even though clustering was only done on locations to 120E
    lons=np.arange(nlon)+cluster_lon_range[0]
    nlat=(cluster_lat_range[1]-cluster_lat_range[0])+1
    lats=np.arange(nlat)+cluster_lat_range[0]
    nrows=2
    ncols=K
    cmap=plt.get_cmap('gnuplot2_r')
    sub_cluster_fig_labels=['a','b','c']
    fig = plt.figure(figsize=(12, 9*nrows/ncols))
    for i in range(2*K):
        if i<K:
            cluster_title='Cl {c:d} '.format(c=i+1)
        else:
            cluster_title='Cl {c1:d}.'.format(c1=sub_cix+1)+sub_cluster_fig_labels[i-K]+' '
        if i==sub_cix:
            ix=np.where(cluster_labels>=K) # find bv indices for the all sub cluster labels
        else:
            ix=np.where(cluster_labels==i) # find bv indices just for this cluster
            cluster_title=cluster_title+cluster_names[i]
        density, this_nbvs=get_track_density(bv_tracks[ix[0]], lons, lats, ndays)
        print(this_nbvs, ' BVs in '+cluster_title)
        title=cluster_title+': {c:d} ({pc:.1f} % BVs)'.format(c=this_nbvs, pc=100*this_nbvs/nbvs)
        con,ax=plot_density(fig, nrows,ncols,i+1, 30*density, density_levels, cmap, title, lons,lats, extend='both')
        this_cluster_centers=cluster_centers[i,:]
        ax.scatter(this_cluster_centers[lon_ix], this_cluster_centers[lat_ix], marker='.', color='cyan',zorder=4)
    cb_pos=[0.1,0.05,0.8,0.02]
    cb = fig.add_axes(cb_pos)
    fig.colorbar(con, cax=cb,orientation = 'horizontal')
    cb.set_xlabel(density_title)
    return plt
    
"""
do_clustering()
This does the actual K means clustering and writes the output to .npy files
It returns the clustered BVs
Inputs:
    bv_tracks850 - the BVs we should potentially cluster (they must meet conditions about being in an area long enough to get nsample track positions)
    outdir - where to store the output
Returns:
    clust_bv_tracks - a subset of bv_tracks850 that were clustered
"""
from sklearn.cluster import KMeans
def do_clustering(bv_tracks850, outdir):
    # To do Kmeans we need nsamples=6 BV locations spread evenly across the time that the BV exists within
    # the usual lon/lat range of a BV (some of them go out of this region
    # so we can only use those bvs which are within the lat/lon range for nsample times
    in_reg=np.asarray([len(bv.exists_in_region(cluster_lon_range, cluster_lat_range)[0]) for bv in bv_tracks850])
    ix=np.where(in_reg>=nsamples)
    print(len(ix[0]), 'are in reg for >=nsamples')

    clust_bv_tracks=bv_tracks850[ix]
    nbvs=len(clust_bv_tracks)
    bv_startdates=np.asarray([bv.track_times[0] for bv in clust_bv_tracks])
    bvids=np.asarray([bv.track_id for bv in clust_bv_tracks])

    # get the data we will need for kmeans
    data=np.zeros((nbvs,nsamples*2)) # to hold lat/lon data passed to kmeans
    for n in range(nbvs):
        tix=clust_bv_tracks[n].exists_in_region(cluster_lon_range, cluster_lat_range)
        lons=clust_bv_tracks[n].vort_data[0,tix[0]]
        lats=clust_bv_tracks[n].vort_data[1,tix[0]]
        this_nt=len(tix[0])
        # find evenly spaced times across lifetime in this region
        rem=this_nt % nsamples
        if rem==0:
            tix=np.arange(nsamples)*int(np.floor(this_nt/nsamples))
        else:
            tix=np.arange(nsamples)*int(np.floor(this_nt/(nsamples-1)))
        if tix[-1]==this_nt:
            tix[-1]=this_nt-1

        data[n,lon_ix]=lons[tix]
        data[n,lat_ix]=lats[tix]

    # now do 1st level of clustering
    model = KMeans(n_clusters=K, init='k-means++', n_init=10, random_state=0) 
    this_kmeans = model.fit(data)
    cluster_centers=np.zeros((K*2,nsamples*2)) # store the centroids for k=3 and subk=3
    cluster_centers[:K,:]=this_kmeans.cluster_centers_
    cluster_labels=this_kmeans.labels_

    # now sub cluster the largest cluster
    unique, counts = np.unique(cluster_labels, return_counts=True)
    sub_cix=np.argmax(counts)
    ix=np.where(cluster_labels==sub_cix)
    data_to_subcluster=data[ix[0],:]

    model = KMeans(n_clusters=K, init='k-means++', n_init=10, random_state=0) 
    this_kmeans = model.fit(data_to_subcluster)
    #overwrite the labels for the subclusters
    cluster_labels[ix[0]]=this_kmeans.labels_+K
    # store these cluster centres in second lot of centres
    cluster_centers[K:,:]=this_kmeans.cluster_centers_

    # now save the output
    np.save(outdir+'kmeans{K:d}_plus_sub.npy'.format(K=K), cluster_labels)
    np.save(outdir+'kmeans{K:d}_plus_sub_centers.npy'.format(K=K), cluster_centers)
    np.save(outdir+'kmeans_bvids.npy', bvids)
    np.save(outdir+'kmeans_bv_startdates.npy', bv_startdates)

    # return the bv_tracks that were used in the clustering
    return clust_bv_tracks