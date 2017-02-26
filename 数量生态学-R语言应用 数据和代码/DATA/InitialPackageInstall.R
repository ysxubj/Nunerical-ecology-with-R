# ============================================================
# Script for first-time users of "Numerical Ecology with R"  #
# by Daniel Borcard, Francois Gillet and Pierre Legendre     #
# ============================================================

# This script installs or provides guidelines to install all 
# the packages necessary to run the code provided in the book, 
# but that do not belong to the standard R distribution (steps 2-3). 

# Steps 1 to 3 must be run only once when installing or upgrading R.
# Step 4 is not mandatory.


# 1. Update installed packages
#    -------------------------
update.packages(checkBuilt=TRUE, ask=FALSE)


# 2. Install packages from the main CRAN site
#    ----------------------------------------

install.packages(c("ade4", "ape", "cluster", "ellipse", "FactoMineR", "FD", 
	"gclus", "labdsv", "MASS", "mvpart", "MVPARTwrap", "RColorBrewer", "spdep", 
	"tripack", "vegan"))
# package "MVPARTwrap" for R 2.15.x is available from CRAN


# 3. Install packages from R-Forge.R-project
#    ---------------------------------------

install.packages("spacemakeR", repos="http://R-Forge.R-project.org")
install.packages("PCNM", repos="http://R-Forge.R-project.org")
install.packages("AEM", repos="http://R-Forge.R-project.org")
# package "packfor" is installed as dependency with "PCNM" for R 2.15.x


# 4. OPTIONAL (for power users): Install all R packages from Environmetrics,
#    a CRAN Task View for the Analysis of Ecological and Environmental Data
#    See http://cran.r-project.org/web/views/Environmetrics.html
#    ----------------------------------------------------------------------

install.packages("ctv")
library(ctv)
install.views("Environmetrics")

# Other potentially useful CRAN Task Views...
install.views("Cluster")
install.views("Multivariate")
install.views("Spatial")
install.views("MachineLearning")
