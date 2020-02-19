# kinase-activity-net


Here is how to make predictions.
    Access predictor files and network  from supplementary tables from Prediction (S1-S5) of signed protein-kinase regulatory circuitss

    The following scripts can be used to make kinase-kinase regulatory relationship prediction
        do-bart-lo-kin.sh
        do-bart-build.sh    to generate bart model
        do-bart-predict.sh   to make predictions
        do-bart-eval.sh     to cross validate model
        do-sign-bart-lo-kin.sh


    For clustering analysis go to src/Netw-clust/
    and type make enrich

    for analysis with independent experimental data go to src/net-validation and type:
    make psites to get differentially phosphorylated phospho sites
    make paths to find shortest paths from inhibited kinase to downregulated phosphosite
    make new_paths to find paths with edges predicted by the in vitro and in-vivo data.
