# faradr

Tools for sequencing data analysis.

Exploration and plotting tools usefull in Next Generation Sequencing data analysis.

## Installation:
The package can be installed with [devtools](https://github.com/hadley/devtools) (by Hadley Wickham)
```
install.package("devtools")
library(devtools)
devtools::install_github("devtools")
```
VLPB github repository is private one, therefore it is neccessery to
authenticate in order to download package. Devtools uses github authentication
tokens. To generate authentication token for R command line follow instructins
from [Github
help](ttps://help.github.com/articles/creating-an-access-token-for-command-line-use)

Install faradr with devtools
```
install_github("UvA-MAD/faradr", auth_user="github_user_id",
        auth_token="github_generated_token")
```
rembember to replace "github_user_id" and "github_generated_token" with your
login credentials.


