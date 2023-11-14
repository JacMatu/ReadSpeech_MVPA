#!/usr/bin/env bash

# for instructions on how to name it according to the lab guide
# see https://github.com/cpp-lln-lab/dataset_maintenance 
# example: 3018_Shire_LOTR_FB_raw
GIN_BASENAME='analyses_MultimodalLanguageMVPA-derivatives'

## CREATE A DATALAD REPO
#datalad root folder must exist 
#cd to inputs folder! 
cd /Volumes/Slim_Reaper/Projects/Language_MVPA/outputs

datalad create --force -c text2git derivatives
cd derivatives


root_dir=${PWD}

preproc_dir=${root_dir}/bidspm-preproc
stats_dir=${root_dir}/bidspm-stats
roi_dir=${root_dir}/bidspm-roi
cosmo_dir=${root_dir}/cosmo-mvpa
# optional

# create the repo on gin for the root directory

#hosted private or cpp-lln-lab organization?
#datalad create-sibling-gin -d . -s origin --access-protocol ssh --private --credential JacMat_gin_token cpp-lln-lab/"${GIN_BASENAME}"
datalad create-sibling-gin -d . -s origin --access-protocol ssh --private --credential JacMat_gin_token "${GIN_BASENAME}"


## create the subdatasets for each analyses step
#PREPROC
datalad create --force -d . -c text2git "${preproc_dir}" 
cd "${preproc_dir}"
datalad create-sibling-gin -d . -s origin --access-protocol ssh --private --credential JacMat_gin_token "${GIN_BASENAME}"-bidspm-preproc

#STATS
cd "${root_dir}"
datalad create --force -d . -c text2git "${stats_dir}" 
cd "${stats_dir}"
datalad create-sibling-gin -d . -s origin --access-protocol ssh --private --credential JacMat_gin_token "${GIN_BASENAME}"-bidspm-stats

#ROI
cd "${root_dir}"
datalad create --force -d . -c text2git "${roi_dir}" 
cd "${roi_dir}"
datalad create-sibling-gin -d . -s origin --access-protocol ssh --private --credential JacMat_gin_token "${GIN_BASENAME}"-bidspm-roi

#MVPA - skipped for now, work in progress! 
cd "${root_dir}"
datalad create --force -d . -c text2git "${cosmo_dir}" 
cd "${cosmo_dir}"
datalad create-sibling-gin -d . -s origin --access-protocol ssh --private --credential JacMat_gin_token "${GIN_BASENAME}"-cosmo-mvpa

#datalad subdatasets --set-property url https://gin.g-node.org/cpp-lln-lab/"${GIN_BASENAME}"-source "${sourcedata}"

datalad subdatasets --set-property url https://gin.g-node.org/JacMatu/"${GIN_BASENAME}"-bidspm-preproc "${bidspm-preproc}"
datalad subdatasets --set-property url https://gin.g-node.org/JacMatu/"${GIN_BASENAME}"-bidspm-stats "${bidspm-stats}"
datalad subdatasets --set-property url https://gin.g-node.org/JacMatu/"${GIN_BASENAME}"-bidspm-roi "${bidspm-roi}"



cd "${root_dir}"

datalad push --to origin -r

echo "############################"
echo "# DATALAD IS READY TO WORK #"
echo "############################"