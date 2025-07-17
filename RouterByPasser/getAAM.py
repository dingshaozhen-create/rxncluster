import subprocess
import os
current_path = os.path.dirname(os.path.abspath(__file__))
import tempfile
import shutil
mediapath = current_path+'/media'
from NeutraliseCharges import NeutraliseCharges as nc
from dszCoder import RetroGenerator
import rdkit.Chem as Chem          # molecule building
from rdkit.Chem import AllChem
from rdkit.Chem import MolFromSmiles, MolToSmiles,MolFromSmarts,rdFMCS,ReplaceSubstructs,ReplaceCore
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit import DataStructs
import json
import re
project_dir = os.path.dirname(current_path)
def createTemp():
	tempath = tempfile.mktemp()
	return tempath.strip().split('/')[-1]   
def mymovefile(srcfile,dstfile):
    if not os.path.isfile(srcfile):
        print ("%s not exist!"%(srcfile))
    else:
        fpath,fname=os.path.split(dstfile)    
        if not os.path.exists(fpath):
            os.makedirs(fpath)                
        shutil.move(srcfile,dstfile)

def getRulesFromSmirks(smirks,timeout=300):
	command = ['java','-jar',current_path+'/aam/RDT1.4.2.jar','-Q','SMI','-q',smirks,'-g','-j','AAM','-f','TEXT']
	try:
		subprocess.run(command, timeout=timeout)
		
	except subprocess.TimeoutExpired:
		return 'AAM is timeout'
	aam_result = os.path.dirname(current_path)+'/ECBLAST_smiles_AAM.txt'

	print (aam_result)
	aam_result = open(aam_result).read().strip().split('//')[1].strip().split("SELECTED AAM MAPPING")[1].lstrip()
	tempath = createTemp()
	savepath = mediapath+'/%s/'%tempath
	mymovefile(os.path.dirname(current_path)+'/ECBLAST_smiles_AAM.png',savepath)
	mymovefile(os.path.dirname(current_path)+'/ECBLAST_smiles_AAM.txt',savepath)
	
	savepath = savepath.replace(project_dir,'')

	rg = RetroGenerator.RetroGenerator()
	print (aam_result)
	smarts = rg.getRetroString(aam_result)

	RadiusDct = {}
	return (aam_result,savepath+'ECBLAST_smiles_AAM.png',smarts)




def reaction_predict(rule,reactants):
	rxn = AllChem.ReactionFromSmarts(rule)
	rxnrs = [Chem.MolFromSmiles(reactant) for reactant in reactants]
	try:
		outcomes = rxn.RunReactants(rxnrs,maxProducts=1000)
	except:
		return []
	result_pro_all = []			
	for outcome in outcomes:
		result_pro = []
		for product in outcome:
			result_pro.append(nc(Chem.MolToSmiles(product))[0])
		result_pro_all.append(tuple(result_pro))	
	return list(set(result_pro_all))

if __name__ == '__main__':
	pass