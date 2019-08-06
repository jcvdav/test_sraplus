# Setup packages

print("############# UPDATING PACKAGES  #############")
## Update existing packages
# update.packages(ask = FALSE, checkBuilt = TRUE, repos = "http://cran.us.r-project.org")
print("############# ALL PACKAGES UPDATED #############")
print("############# INSTALLING DEVTOOLS  #############")
## Check for devtools
if (!require(devtools)){
  install.packages(pkgs = "devtools", repos = "http://cran.us.r-project.org")
}
## Install sraplus
remotes::install_github("danovando/sraplus", ref = 'v2.0')

# Run the test
startR::render_doc(file = "test.Rmd", rmd_dir = ".", output_dir = ".")
