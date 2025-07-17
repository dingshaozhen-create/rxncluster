from django.db import models

# Create your models here.

class Molecule(models.Model): 
	dbid = models.CharField(max_length=50)
	smi = models.CharField(max_length=2000)
	smi_nc = models.TextField()
	dblabel = models.TextField()

class molecule_id_name(models.Model):
	dbid = models.CharField(max_length = 50)
	name = models.CharField(max_length = 1000)
	namelength = models.IntegerField()

class bbs_rxn_merge(models.Model):
	bbsRxnID = models.IntegerField()
	firstcompdsmi = models.CharField(max_length=2000)
	endcompdsmi = models.CharField(max_length=2000)
	smirks = models.TextField()
	mapped_smirks = models.CharField(max_length = 3000)
	smarts = models.CharField(max_length = 3000)
	simplesmarts = models.CharField(max_length=3000)
	bbsretroid = models.IntegerField()

class bbs_single_step(models.Model):
	rxn = models.CharField(max_length = 50)
	reactant = models.CharField(max_length=50)
	product = models.CharField(max_length=50)
	reactant_smi = models.CharField(max_length=2000)
	product_smi=models.CharField(max_length=2000)
	smirks=models.CharField(max_length=2000)
	#mapped_smirks = models.CharField(max_length=3000)
	bbsRxnID = models.IntegerField()

class rhea_kegg_ec(models.Model):
	rxnid = models.CharField(max_length=20)
	direction = models.CharField(max_length=20)
	master_id = models.CharField(max_length=20)
	ec = models.CharField(max_length=50)


class pathsearch_combo_steps_smarts(models.Model):
	step1 = models.CharField(max_length=50)
	step2 = models.CharField(max_length=50)
	step3 = models.CharField(max_length=50)
	step4 = models.CharField(max_length=50)

	reactant_1 = models.CharField(max_length=50)
	product_1 = models.CharField(max_length=50)
	reactant_2 = models.CharField(max_length=50)
	product_2 = models.CharField(max_length=50)
	reactant_3 = models.CharField(max_length=50)
	product_3 = models.CharField(max_length=50)
	reactant_4 = models.CharField(max_length=50)
	product_4 = models.CharField(max_length=50)

	steplength = models.IntegerField()
	bbsRxnID = models.IntegerField()

	step1_smirks = models.TextField()
	step2_smirks = models.TextField()
	step3_smirks = models.TextField()
	step4_smirks = models.TextField()

	smarts1 = models.CharField(max_length=3000)
	bbsRxnID_1 = models.IntegerField()

	smarts2 = models.CharField(max_length=3000)
	bbsRxnID_2 = models.IntegerField()

	smarts3 = models.CharField(max_length=3000)
	bbsRxnID_3 = models.IntegerField()

	smarts4 = models.CharField(max_length=3000)
	bbsRxnID_4 = models.IntegerField()

	combo_simple_smarts_count = models.IntegerField()
	combo_simple_smarts = models.CharField(max_length=3000)

	idx = models.IntegerField()

	






# class rxnCombo(models.Model):
# 	rxnCombo = models.CharField(max_length = 2000)
# 	step1 = models.CharField(max_length = 50)
# 	step2 = models.CharField(max_length = 50)
# 	step3 = models.CharField(max_length = 50)
# 	step4 = models.CharField(max_length = 50)
# 	reactant_1 = models.CharField(max_length = 50)
# 	product_1 = models.CharField(max_length=50)
	
# 	reactant_2 = models.CharField(max_length = 50)
# 	product_2 = models.CharField(max_length=50)
	
# 	reactant_3 = models.CharField(max_length = 50)
# 	product_3 = models.CharField(max_length=50)

# 	reactant_4 = models.CharField(max_length = 50)
# 	product_4 = models.CharField(max_length=50)

# 	firstcompd = models.CharField(max_length = 50)
# 	endcompd = models.CharField(max_length=50)

# 	steplength = models.IntegerField()

# 	firstcompdsmi =models.CharField(max_length = 2000)
# 	endcompdsmi = models.CharField(max_length=2000)

# 	bbsRxnID = models.IntegerField()

# class retroRules(Models.model):
# 	bbsRxnID = models.IntegerField()
# 	smarts = models.CharField(max_length=2000)
# 	smirks = models.CharField(max_length = 2000)

# class mets(Models.model):
# 	dbid = models.CharField(max_length = 50)
# 	smi = models.CharField(max_length = 2000)
# 	smi_nc = models.CharField(max_length = 2000)
# 	dblabel = models.CharField(max_length)







