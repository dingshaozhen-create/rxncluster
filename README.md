## Overview

This repository contains the source code for RxnCluster, a web-based tool designed to explore reaction clusters leading to target molecules by digitalizing typical biosynthetic patterns(website: http://design.rxnfinder.org/rxncluster/).


## Technical Architecture
RxnCluster is built using:
Backend: Django(Python web framework,version=3.2.23)

Database: PostgreSQL (for structured data storage)

Frontend: HTML, CSS, JavaScript with Django templates

Chemistry Toolkit: RDKit

Visualization: Custom JavaScript components for reaction cluster display

## Steps
### 1. Clone the repository:
   ```
   git clone https://github.com/dingshaozhen-create/rxncluster.git
   cd rxncluster
   ```

### 2. Set up the environment
We strongly recommand to install packages in [conda](https://docs.conda.io/projects/miniconda/en/latest/miniconda-install.html) environment.

After install conda, create an environmnet
```
conda create -n rxncluster python==3.7
```
activate the environment and install the packages

```
conda activate rxncluster
cd /dir/to/this/repo # !cd the dir you git clone or download this repository
pip install -r requirements.txt
```

### 3. Load the data:

All data utilized in this project is stored in /path/to/rxncluste/DataExcel/. These Excel files need to be imported into the PostgreSQL database. We recommend using Navicat (https://www.navicat.com.cn/) as the database management tool to import the data into PostgreSQL.

### 4. Set up the PostgreSQL database:
   - Create a database named `rxncluster`.
   - Update the database configuration in `settings.py`:
     
     ```python
     DATABASES = {
         'default': {
             'ENGINE': 'django.db.backends.postgresql_psycopg2',
             'NAME': 'rxncluster',
             'USER': 'your_username',
             'PASSWORD': 'your_password',
             'HOST': 'localhost',
             'PORT': '5432',
         }
     }
     ```


### 5. How to Run

1. Start the Django development server:
   ```
   conda activate rxncluster
   cd /path/to/rxncluster/
   python manage.py runserver
   
   ```
2. Open a web browser and navigate to:
   ```
   http://127.0.0.1:8000/rxncluster/
   ```

## Using RxnCluster by WebServer
Please see http://design.rxnfinder.org/rxncluster/faq/
RxnCluster provides API for users who are not good at computer programming.

## Using RxnCluster by API

  To get complete result of reaction clusters
   ```
   http://design.rxnfinder.org/rxncluster/API_result_request/precursor_smiles/target_smiles
   ```
   
   Example: 
   ```
   http://design.rxnfinder.org/rxncluster/API_result_request/NC(Cc1cccc2ccccc12)C(=O)O/O=C(O)Cc1cccc2ccccc12
   ```

   To get the reaction rule 
   
   ```
   http://design.rxnfinder.org/rxncluster/API_rule_extract/{reactant}/{product}
   ```
  Example:
  ```
  http://design.rxnfinder.org/rxncluster/API_rule_extract/NC(Cc1cccc2ccccc12)C(=O)O/O=C(O)Cc1cccc2ccccc12
   ```

## License
This project is licensed under the MIT License.

## Citation
If you use RxnCluster in your research, please cite the original paper:

```
Ding et al., "RxnCluster: A web-based tool for exploring reaction clusters leading to target molecules by digitalizing typical biosynthetic patterns", ACS Synthetic Biology, 2025.
```

## Contact
For questions or feedback, please contact:
- Shaozhen Ding [23113176@whpu.edu.cn]
- Project repository: [https://github.com/dingshaozhen-create/rxncluster/]