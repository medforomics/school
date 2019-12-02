#!/bin/bash


module load azure/2.0.72 
#export HTTP_PROXY="http://proxy.swmed.edu:3128"
#export HTTPS_PROXY="http://proxy.swmed.edu:3128"
#export NO_PROXY="127.0.0.1,localhost"
#export http_proxy="http://proxy.swmed.edu:3128"
#export https_proxy="https://proxy.swmed.edu:3128"
export AZURE_STORAGE_ACCOUNT=swazrstrseq
export AZURE_STORAGE_KEY=SxfAj0wkNDyQVKcP5ChoTDEd7J8bZm2zAqh0vNE2YRAQxFGrTLc1wvlWv0IYrS1p6thpBCvBLGHbVQmiu1/XqQ==

#list containers
az storage container list

#list files in container
az storage blob list -c container_name

# upload files from a directory to a blob container
az storage blob upload-batch -d container -s directory
