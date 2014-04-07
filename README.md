# faradr

Tools for sequencing data analysis.

Exploration and plotting tools usefull in Next Generation Sequencing data analysis.

## Installation:
The package can be installed with [devtools](https://github.com/hadley/devtools) (by Hadley Wickham)
```R
install.packages("devtools")
library(devtools)

# if using latest R version (>= 3.0.3)
# update devtools to the latest vesion
devtools::install_github("devtools")

# check version of devtools installed
packageVersion("devtools")
```
VLPB github repository is private one, therefore it is neccessery to
authenticate in order to download package. 

**devtools >= 1.5** uses github authentication tokens. To generate authentication token for R command line follow instructions
from [Github
help](https://help.github.com/articles/creating-an-access-token-for-command-line-use)
Install faradr with devtools
```R
install_github("UvA-MAD/faradr", auth_user="github_user_id",
        auth_token="github_generated_token")
```
rembember to replace "github_user_id" and "github_generated_token" with your
login credentials.

**devtools < 1.5** uses github password to authenticate. 
```R
install_github("UvA-MAD/faradr", auth_user="github_user_id",
        password="github_account_password")
```
rembember to replace "github_user_id" and "github_account_password" with your
login credentials.