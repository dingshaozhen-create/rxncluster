import dszCode
from django.shortcuts import render
from rdkit import Chem
from django.http import HttpResponse,JsonResponse
from NeutraliseCharges import NeutraliseCharges as nc
import json
from .getAAM import getRulesFromSmirks
import rdkit
import psycopg2
import re
import tempfile
import os
current_path = os.path.dirname(os.path.abspath(__file__))
import shutil
from rdkit.Chem import AllChem, Draw, MolFromSmiles, MolFromSmarts
from rdkit.Chem import rdMolHash
molhashf = rdMolHash.HashFunction.names
from io import StringIO,BytesIO
from PIL import Image, ImageDraw, ImageFont
import textwrap
import os.path
import urllib.parse
from rdkit.Chem.Draw import rdMolDraw2D
from dszCoder import RetroGenerator,calBulkSimilarity
import pandas as pd
import math
from django.db.models import Q
mediapath = current_path+'/media'
from django.core.paginator import Paginator,PageNotAnInteger,EmptyPage
from django.template.loader import render_to_string
from RouterByPasser.models import Molecule,bbs_rxn_merge,bbs_single_step,bbs_single_step,rhea_kegg_ec,pathsearch_combo_steps_smarts,molecule_id_name

#bbsRxnIDsRemove = [133215]

def createTemp():
	tempath = tempfile.mktemp()
	return tempath.strip().split('/')[-1]  
def rxnCluster(request):
	return render(request,"rxnCluster.html")

def convert_marksmiles_to_image_stream(smiles,smarts,size=(300,300)):
	mol = Chem.MolFromSmiles(smiles)
	patt = Chem.MolFromSmarts(smarts)
	hit_at = mol.GetSubstructMatch(patt)
	hit_bond = []
	for bond in patt.GetBonds():
		aid1 = hit_at[bond.GetBeginAtomIdx()]
		aid2 = hit_at[bond.GetEndAtomIdx()]
		hit_bond.append(mol.GetBondBetweenAtoms(aid1, aid2).GetIdx())

	img = Draw.MolToImage(mol, size=size,highlightAtoms=list(hit_at),highlightBonds=hit_bond)
	buf = BytesIO()
	img.save(buf, "PNG")
	image_stream = buf.getvalue()
	return image_stream
def convert_smiles_to_image_stream(smiles, size=(300, 300)):
	mol = AllChem.MolFromSmiles(str(smiles))
	img = Draw.MolToImage(mol, size=size)
	buf = BytesIO()
	img.save(buf, "PNG")
	image_stream = buf.getvalue()
	return image_stream

def convert_retro_to_image_stream(smiles, size=(300, 600)):
	rxn = AllChem.ReactionFromSmarts(str(smiles))
	img = Draw.ReactionToImage(rxn, tuple(size))
	buf = BytesIO()
	img.save(buf, "PNG")
	image_stream = buf.getvalue()
	return image_stream

def png_response(data):
	response = HttpResponse(data, content_type="image/png")
	response['Last-Modified'] = 'Mon, 27 Apr 2015 02:05:03 GMT'
	response['Cache-Control'] = 'max-age=31536000'
	return response	
def marksmiles2img(request,smi,smarts,width,height):
	width = int(width)
	height = int(height)
	data = convert_marksmiles_to_image_stream(smi, smarts,size=(width,height))
	return png_response(data)
def smi2img(request,smi,width,height):
	width = int(width)
	height = int(height)
	data = convert_smiles_to_image_stream(smi, size=(width,height))
	return png_response(data)

def retro2img(request,rules,width,height):
	width = int(width)
	height = int(height)
	data = convert_retro_to_image_stream(rules, size=(width,height))
	return png_response(data)

def smiles_to_mol(request):
	smiles = request.GET.get('smiles')
	try:
		m = Chem.MolFromSmiles(smiles)
		mol = Chem.MolToMolBlock(m)
	except:
		mol='error'
	return HttpResponse(json.dumps({"mol":mol}),content_type='application/json')

def sortedSingleCluster(single_cluster):

	df = pd.DataFrame(single_cluster)
	sorted_df = df.sort_values(by='smiScore', ascending=False)
	result = sorted_df.to_dict("list")
	return result

def getExactMassAndFormula(smi):
	from rdkit.Chem import Descriptors
	mol = Chem.MolFromSmiles(smi)
	exact_mass = round(Descriptors.ExactMolWt(mol),3)
	formula = rdMolHash.MolHash(mol, rdkit.Chem.rdMolHash.HashFunction.MolFormula)
	return (exact_mass,formula)

def getNameFromDBid(dbid):
	dbid = dbid.replace("KEGG:","").replace("ChEBI:","").replace("BiGG:","")
	namesobjs = molecule_id_name.objects.filter(dbid=dbid).order_by("namelength")
	objTemp = ''
	if len(namesobjs) == 0:
		return ""
	for obj in namesobjs:
		objTemp = obj
		if obj.name != None and len(obj.name)>=3:
			return obj.name
	return objTemp.name


def getResultPrecursorTarget(precursor_struc,target_struc):
	precursor_struc = nc(precursor_struc.replace("\\","").replace("/",''))[0]
	target_struc = nc(target_struc.replace("\\","").replace("/",''))[0]
	rxnSmirks = precursor_struc+">>"+target_struc
	smarts_result = getRulesFromSmirks(rxnSmirks)
	aamResult = smarts_result[0]
	png_path = smarts_result[1]
	retroString = smarts_result[2]
	result_total = getKnownRouteFromRetro(retroString,precursor_struc,target_struc)
	result_total_new = []
	for i in result_total:
		i["Score"] = round(math.log(i["combo_simple_smarts_count"])+(10-sum(i["SaScore"])/i["steplength"])+(4-i["steplength"]),2) 
		result_total_new.append(i)

	result_total = result_total_new

	result_total1 = sorted([i for i in result_total if i["steplength"] == 1],key=lambda x:x["Score"],reverse=True)
	
	
	result_total2 = sorted([i for i in result_total if i["steplength"] == 2],key=lambda x:x["Score"],reverse=True)
	result_total3 = sorted([i for i in result_total if i["steplength"] == 3],key=lambda x:x["Score"],reverse=True)
	result_total4 = sorted([i for i in result_total if i["steplength"] == 4],key=lambda x:x["Score"],reverse=True)

	result_total1 = merge_result_total_n(result_total1)
	result_total2 = merge_result_total_n(result_total2)
	result_total3 = merge_result_total_n(result_total3)
	result_total4 = merge_result_total_n(result_total4)

	for each_one in result_total1:
		each_one["single_cluster"] = sortedSingleCluster(each_one["single_cluster"])
	for each_one in result_total2:
		each_one["single_cluster"] = sortedSingleCluster(each_one["single_cluster"])
	for each_one in result_total3:
		each_one["single_cluster"] = sortedSingleCluster(each_one["single_cluster"])
	for each_one in result_total4:
		each_one["single_cluster"] = sortedSingleCluster(each_one["single_cluster"])

	tempath = createTemp()
	savepath = mediapath+'/%s.json'%tempath
	result = {'step1':result_total1,'step2':result_total2,'step3':result_total3,'step4':result_total4}
	fwrite = open(savepath,'w')
	json.dump(result,fwrite)
	fwrite.close()

	new_dict = {"aamResult":aamResult,"png_path":png_path,"retroString":retroString,"result_total1":result_total1,"result_total2":result_total2,"result_total3":result_total3,"result_total4":result_total4,"tempath":tempath}

	return new_dict

def search(request):
	precursor_struc = request.GET.get("precursor_struc")
	target_struc = request.GET.get("target_struc")
	if precursor_struc == "":
		similarCompd = calBulkSimilarity.calculateSimilarities(current_path+"/static/endcompdsmis.txt",target_struc,10,0)
		score = [round(i[1],3) for i in similarCompd]
		compdsmis = [i[0] for i in similarCompd]
		dblinks = [getdbidFromSmi(smi) for smi in compdsmis]
		names = [getNameFromDBid(dbid) for dbid in dblinks]
		exactMassAndFormula = [getExactMassAndFormula(smi) for smi in compdsmis]
		exactMass = [i[0] for i in exactMassAndFormula]
		formula = [i[1] for i in exactMassAndFormula]

		return render(request,"rxncluster_similarity.html",{"score":score,"compdsmis":compdsmis,"target_struc":target_struc,"dblinks":dblinks,"exactMass":exactMass,"formula":formula,"names":names})
	
	new_dict = getResultPrecursorTarget(precursor_struc,target_struc)
	
	return render(request,"rxnCluster_result.html",new_dict)


def getSimpleSmarts(retroString):
	rg = RetroGenerator.RetroGenerator()
	retroStringList = retroString.split("$$")
	return "$$".join([rg.getSimpleSmarts(retroString) for retroString in retroStringList])

def result(request,tempfilename):
	result_total = json.load(open(current_path+'/media/tempfile/%s.json'%tempfile_name))
	result_total1 = [i for i in result_total if i["steplength"] == 1]
	result_total2 = [i for i in result_total if i["steplength"] == 2]
	result_total3 = [i for i in result_total if i["steplength"] == 3]
	result_total4 = [i for i in result_total if i["steplength"] == 4]

	
	return render(request,"rxnCluster_result.html",{"aamResult":aamResult,"png_path":png_path,"retroString":retroString,"result_total1":result_total1,"result_total2":result_total2,"result_total3":result_total3,"result_total4":result_total4,"tempath":tempath})

def getKnownRouteFrombbsRxnIDs(bbsRxnIDs,precursor_struc,target_struc):
	result_sql3 = []
	for bbsRxnID in bbsRxnIDs:
		result_sql3.extend(pathsearch_combo_steps_smarts.objects.filter(bbsRxnID=bbsRxnID))

	result_total = []
	
	noresultSimpleSmarts = []


	for each_one in result_sql3:
		combo_simple_smarts_count = each_one.combo_simple_smarts_count
		if combo_simple_smarts_count == 0:
			continue
		combo_simple_smarts = each_one.combo_simple_smarts
		bbsRxnID = each_one.bbsRxnID
		step1 = each_one.step1
		step2 = each_one.step2
		step3 = each_one.step3
		step4 = each_one.step4
		reactant1 = each_one.reactant_1
		product1 = each_one.product_1
		reactant2 = each_one.reactant_2
		product2 = each_one.product_2
		reactant3 = each_one.reactant_3
		product3 = each_one.product_3
		reactant4 = each_one.reactant_4
		product4 = each_one.product_4
		steplength = each_one.steplength
		smarts1 = each_one.smarts1
		smarts2 = each_one.smarts2
		smarts3 = each_one.smarts3
		smarts4 = each_one.smarts4

		smarts_reactant_1 = ''
		smarts_product_1 = ''

		smarts_reactant_2 = ''
		smarts_product_2 = ''

		smarts_reactant_3 = ''
		smarts_product_3 = ''

		smarts_reactant_4 = ''
		smarts_product_4 = ''


		if smarts1 != None:
			smarts_reactant_1 = smarts1.split(">>")[0]
			smarts_product_1 = smarts1.split(">>")[1]

		if smarts2 != None:
			smarts_reactant_2 = smarts2.split(">>")[0]
			smarts_product_2 = smarts2.split(">>")[1]

		if smarts3 != None:
			smarts_reactant_3 = smarts3.split(">>")[0]
			smarts_product_3 = smarts3.split(">>")[1]

		if smarts4 != None:
			smarts_reactant_4 = smarts4.split(">>")[0]
			smarts_product_4 = smarts4.split(">>")[1]

		rxns = ';'.join([i for i in [step1,step2,step3,step4] if i != 'norxn']);
		compds = ''
		if steplength == 1:
			compds = [reactant1,product1]
			smartsList = [smarts_reactant_1,smarts_product_1]
		elif steplength == 2:
			compds = [reactant1,product1,product2]
			smartsList = [smarts_reactant_1,smarts_product_1,smarts_product_2]
		elif steplength == 3:
			compds = [reactant1,product1,product2,product3]
			smartsList = [smarts_reactant_1,smarts_product_1,smarts_product_2,smarts_product_3]
		elif steplength == 4:
			compds = [reactant1,product1,product2,product3,product4]
			smartsList = [smarts_reactant_1,smarts_product_1,smarts_product_2,smarts_product_3,smarts_product_4]


		
		comboretro = "$$".join([i for i in [smarts1,smarts2,smarts3,smarts4] if i !=None])
		comboretroUrl ="$$".join([urllib.parse.quote(i) for i in [smarts1,smarts2,smarts3,smarts4] if i != None])

		simpleComboretro = getSimpleSmarts(comboretro)

		if simpleComboretro in noresultSimpleSmarts:
			continue
		
		each_detail = {}
		each_detail["step1"] = step1
		each_detail["step2"] = step2
		each_detail["step3"] = step3
		each_detail["step4"] = step4
		each_detail["reactant1"] = reactant1
		each_detail["product1"] = product1

		each_detail["reactant2"] = reactant2
		each_detail["product2"] = product2

		each_detail["reactant3"] = reactant3
		each_detail["product3"] = product3

		each_detail["reactant4"] = reactant4
		each_detail["product4"] = product4
		each_detail["rxns"] = rxns
		#each_detail["compds"] = compds
		#each_detail["compds"] = [getCompdSmi(i)[1] for i in compds]

		refProductsList =[Molecule.objects.filter(dbid=i).first() for i in [product1,product2,product3,product4] if i!='nomet']
		refProductsList = [i.smi_nc for i in refProductsList if i != None]
		product_pre_list = getPredictProducts(comboretro,precursor_struc,target_struc,refProductsList)
		
		SaScore = [round(dszCode.getSaScore(i),3) for i in product_pre_list]

		if len(product_pre_list)==0:
			noresultSimpleSmarts.append(simpleComboretro)
			continue
		
		if len([i for i in result_total if i["simplesmarts"] == simpleComboretro]) !=0:
			obj = [i for i in result_total if i["simplesmarts"] == simpleComboretro][0]
			if comboretro not in obj["comboretro"]:
				obj["comboretro"].append(comboretro)
				obj["comboretroUrl"].append(comboretroUrl)
			obj["detail"].append(each_detail)
			compdsmisTemp = [getCompdSmi(i) for i in compds]
			
			if tuple([urllib.parse.quote(i) for i in compdsmisTemp]) not in obj["single_cluster"]["compdsmis"]:
				obj["single_cluster"]["compdsmis"].append(tuple([urllib.parse.quote(i) for i in compdsmisTemp]))
				dbidsTemp = [getdbidFromSmi(i) for i in compdsmisTemp]
				obj["single_cluster"]["compds"].append(dbidsTemp)
				obj["single_cluster"]["smartsList"].append([urllib.parse.quote(i) for i in smartsList])
				obj["single_cluster"]["bbsRxnID"].append(bbsRxnID)
				obj["single_cluster"]["compdnames"].append([getNameFromDBid(i) for i in dbidsTemp])

				obj["single_cluster"]["smiScore"].append(round(dszCode.calculate_similarity_customized(compdsmisTemp[-1],target_struc),2))



		else:
			first_one = {}
			first_one["comboretro"] = [comboretro]
			first_one["comboretroUrl"] =[comboretroUrl]
			first_one["simplesmarts"] = simpleComboretro
			first_one["detail"] = [each_detail]
			first_one['combo_simple_smarts'] = combo_simple_smarts
			first_one["combo_simple_smarts_count"] = combo_simple_smarts_count
			first_one["single_cluster"]  = {}
			compdsmis = [getCompdSmi(i) for i in compds]
			first_one["single_cluster"]["smiScore"] = [round(dszCode.calculate_similarity_customized(compdsmis[-1],target_struc),2)]

			
			
			first_one["single_cluster"]["compdsmis"] = [tuple([urllib.parse.quote(i) for i in compdsmis])]
			first_one["single_cluster"]["smartsList"] = [[urllib.parse.quote(i) for i in smartsList]]

			first_one["single_cluster"]["compds"] = [[getdbidFromSmi(i) for i in compdsmis]]
			first_one["single_cluster"]["compdnames"] =[[getNameFromDBid(i) for i in compds]]
			first_one["product_pre_list"] =[urllib.parse.quote(i) for i in product_pre_list]
			first_one["SaScore"] = SaScore
			first_one["steplength"] = steplength
			first_one["single_cluster"]["bbsRxnID"] = [bbsRxnID]
			result_total.append(first_one)
	for each_one in result_total:
	 	each_one["single_cluster"]["compdsmis"] = [list(i) for i in each_one["single_cluster"]["compdsmis"]]
	return result_total

def getKnownRouteFromRetro(retroString,precursor_struc,target_struc):
	simpleSmarts = getSimpleSmarts(retroString)
	bbsRxnIDs = list(set([i.bbsRxnID for i in bbs_rxn_merge.objects.filter(simplesmarts=simpleSmarts)]))
	#bbsRxnIDs = list(set(bbsRxnIDs)-set(bbsRxnIDsRemove))
	result_total = getKnownRouteFrombbsRxnIDs(bbsRxnIDs,precursor_struc,target_struc)
			
	return result_total


def getPredictProducts(comboretro,precursor_struc,target_struc,refproductsList):
	comboretroList = comboretro.split("$$")
	index = 0
	product_pre_list = [precursor_struc]
	simiScore = round(dszCode.calculate_similarity_customized(target_struc,refproductsList[len(comboretroList)-1]),2)
	simiCompdFroCompare = target_struc
	for index,retro in enumerate(comboretroList):
		reactant_smarts = retro.split(">>")[0]
		product_smarts = retro.split(">>")[1]
		patt = Chem.MolFromSmarts(reactant_smarts)
		mol = Chem.MolFromSmiles(precursor_struc)
		if mol ==None:
			return []
		if not mol.HasSubstructMatch(patt):
			return []
		retro = "("+reactant_smarts +")"+">>"+product_smarts
		rxn = AllChem.ReactionFromSmarts(retro)
		reactants = (AllChem.MolFromSmiles(precursor_struc),)
		products = rxn.RunReactants(reactants)
		if simiScore>=0.7:
			simiCompdFroCompare = refproductsList[index]
		product_score_dict = {nc(Chem.MolToSmiles(i[0]))[0]:round(dszCode.calculate_similarity_customized(nc(Chem.MolToSmiles(i[0]))[0],simiCompdFroCompare),2) for i in products}
		precursor_struc = max(product_score_dict,key=lambda x:product_score_dict[x])
		#precursor_struc = nc(Chem.MolToSmiles(products[0][0]))[0]
		product_pre_list.append(precursor_struc)
	if product_pre_list[-1] != target_struc:
		return []
	return product_pre_list

def showDetail(request):

	tempath = request.GET.get("tempath")

	steplength = request.GET.get("steplength")
	index = int(request.GET.get('index'))
	result = json.load(open(mediapath+"/%s.json"%tempath))['step'+steplength][index]
	combo_simple_smarts_count = result["combo_simple_smarts_count"]
	SaScores = result["SaScore"]
	steplength = result["steplength"]
	Score = math.log(combo_simple_smarts_count)+(10-sum(SaScores)/steplength)+(4-steplength)
	result["Score"] = round(Score,2)
	
	return render(request,"rxncluster_single_cluster.html",{'result':result})


def getCompdSmi(dbid):
	compdsObjs = Molecule.objects.filter(dbid=dbid)
	smi_ncs = [obj.smi_nc for obj in compdsObjs]
	smi = smi_ncs[0]
	return smi

def getdbidFromSmi(smi):
	dict_order = {"kegg":1,"chebi":2,"bigg":3}
	molecules= sorted([(molObj.dbid,molObj.dblabel) for molObj in Molecule.objects.filter(smi_nc=smi)],key=lambda x:dict_order[x[1]])
	dblabel = molecules[0][1]
	dbid = molecules[0][0]
	if dblabel == 'kegg':
		dblabel ='KEGG'
	elif dblabel =='chebi':
		dblabel = 'ChEBI'
	elif dblabel=='bigg':
		dblabel = 'BiGG'


	return dblabel+":"+dbid

def getRefRxnFromSmirks(request):
	#cur = conn_pg.cursor()
	reactant = request.GET.get("reactant")
	product = request.GET.get("product")
	reactant = urllib.parse.unquote(reactant)
	product = urllib.parse.unquote(product)

	reactant_id = request.GET.get("reactant_id").replace("KEGG:","").replace("ChEBI:","").replace("BiGG:","")
	product_id = request.GET.get("product_id").replace("KEGG:","").replace("ChEBI:","").replace("BiGG:","")

	
	reactionList = list(set([i.rxn for i in bbs_single_step.objects.filter(reactant_smi=reactant,product_smi=product)]))

	if len(reactionList) == 0:

		reactant_smi = Molecule.objects.filter(dbid=reactant_id).first().smi_nc
		product_smi = Molecule.objects.filter(dbid=product_id).first().smi_nc
		reactionList = list(set([i.rxn for i in bbs_single_step.objects.filter(reactant_smi=reactant_smi,product_smi=product_smi)]))
	
	ecList = getECByRxnID([i.replace("add_rxn_kegg_","").replace("add_rxn_","") for i in reactionList])
	return HttpResponse(json.dumps([reactionList,ecList]),content_type='applications/json')

def getECByRxnID(rxnids):
	return list(set([i.ec for i in rhea_kegg_ec.objects.filter(rxnid__in=rxnids)]))
def browse(request):
	after_range_num = 5  # 当前页前显示5页
	befor_range_num = 4  # 当前页后显示4页

	obj_list = [{"name":"xiaoMing"*100,"age":"18"*1000}]*2000000
	paginator = Paginator(obj_list,10)
	page = int(request.GET.get('page',1))
	
	if page <=5:
		after_range_num = page
		befor_range_num = 10-page
	if page >= after_range_num:
		page_range = paginator.page_range[page - after_range_num:page + befor_range_num]
	else:
		page_range = paginator.page_range[0:int(page) + befor_range_num]
	try:
		objects = paginator.get_page(page)
	except PageNotAnInteger:
	 	objects = paginator.page(1) 
	except EmptyPage:
	 	objects = paginator.page(paginator.num_pages)
	if request.is_ajax():
		rendered = render_to_string('rxnCluster_browse.html', {'objects':objects,"page_range":page_range})
		#return HttpResponse(json.dumps({"objects":[obj for obj in objects],"page":objects.number}),content_type='application/json')
		return JsonResponse({
			'html': rendered,
			'page': objects.number,
		})

	return render(request, 'rxnCluster_browse.html', {'objects': objects,"page_range":page_range})

"The reaction-template clutster that have been used in pathway for"

def statistics(request):
	return render(request,"bbs_statistics.html")
def faq(request):
	return render(request,"rxncluster_faq.html")
def getPredictPrecursorFromTarget(comboretro,target_struc):
	'''从一个targetStruc逆向退出precursor_struc,comboretro为组合模板
	'''
	precursors_pre_list = []
	reactant_retro,product_retro = comboretro.strip().split(">>")
	mol = AllChem.MolFromSmiles(target_struc)
	patt = AllChem.MolFromSmarts(product_retro)
	if not mol.HasSubstructMatch(patt):
			return ""
	retro = "("+product_retro+")"+">>"+"("+reactant_retro+")"
	rxn = AllChem.ReactionFromSmarts(retro)
	products = (AllChem.MolFromSmiles(target_struc),)
	precursors = rxn.RunReactants(products)
	precursor_struc = nc(Chem.MolToSmiles(precursors[0][0]))[0]
	return precursor_struc

def markCompd(dbid):
	import re
	pattern1 = r'^\d+$'
	if re.search(pattern1,dbid) is not None:
		return 'ChEBI:'+dbid
	elif re.search(r'^\w\d{5}$',dbid) is not None:
		return 'KEGG:'+dbid
	else:
		return 'BiGG:'+dbid
def getExtendResult_total(result_total_forone,target_struc):

	objstemp = pathsearch_combo_steps_smarts.objects.filter(combo_simple_smarts=result_total_forone["combo_simple_smarts"])
	single_cluster = result_total_forone["single_cluster"]
	compds = single_cluster["compds"]
	compdsmis = single_cluster["compdsmis"]
	smartsList = single_cluster["smartsList"]
	smiScore = single_cluster["smiScore"]
	compdnames = single_cluster["compdnames"]
	idxs = []
	bbsRxnID = single_cluster["bbsRxnID"]
	Score = round(math.log(result_total_forone["combo_simple_smarts_count"])+(10-sum(result_total_forone["SaScore"])/result_total_forone["steplength"])+(4-result_total_forone["steplength"]),2)
	result_total_forone["Score"] = Score
	for each_one in objstemp:
		if each_one.bbsRxnID in bbsRxnID:
			continue
		if each_one.idx != 1:
			continue
		steplength = each_one.steplength
		bbsRxnID.append(each_one.bbsRxnID)

		reactant_1 = markCompd(each_one.reactant_1)
		product_1 = markCompd(each_one.product_1)
		product_2 = markCompd(each_one.product_2)
		product_3 = markCompd(each_one.product_3)
		product_4 = markCompd(each_one.product_4)
		reactant_1_smi,product_1_smi = each_one.step1_smirks.split(">>")
		reactant_1_smarts,product_1_smarts = each_one.smarts1.split(">>")
		if each_one.step2_smirks != None and ">>" in each_one.step2_smirks:
			product_2_smi = each_one.step2_smirks.split(">>")[-1]
			product_2_smarts = each_one.smarts2.split(">>")[-1]
		if each_one.step3_smirks != None and ">>" in each_one.step3_smirks:
			product_3_smi = each_one.step3_smirks.split(">>")[-1]
			product_3_smarts = each_one.smarts3.split(">>")[-1]
		if each_one.step4_smirks != None and ">>" in each_one.step4_smirks:
			product_4_smi = each_one.step4_smirks.split(">>")[-1]
			product_4_smarts = each_one.smarts4.split(">>")[-1]
		if steplength == 1:
			dbidTemp = [reactant_1,product_1]
			compds.append(dbidTemp)
			compdnames.append([getNameFromDBid(i) for i in dbidTemp])
			compdsmis.append([reactant_1_smi,product_1_smi])
			smartsList.append([reactant_1_smarts,product_1_smarts])
			smiScore.append(round(dszCode.calculate_similarity_customized(target_struc,product_1_smi),2))

		elif steplength	== 2:
			dbidTemp = [reactant_1, product_1,product_2]
			compds.append(dbidTemp)
			compdnames.append([getNameFromDBid(i) for i in dbidTemp])

			compdsmis.append([reactant_1_smi, product_1_smi,product_2_smi])
			smartsList.append([reactant_1_smarts, product_1_smarts,product_2_smarts])
			smiScore.append(round(dszCode.calculate_similarity_customized(target_struc, product_2_smi),2))

		elif steplength == 3:
			dbidTemp = [reactant_1, product_1, product_2,product_3]
			compds.append(dbidTemp)
			compdnames.append([getNameFromDBid(i) for i in dbidTemp])

			compdsmis.append([reactant_1_smi, product_1_smi, product_2_smi,product_3_smi])
			smartsList.append([reactant_1_smarts, product_1_smarts, product_2_smarts,product_3_smarts])
			smiScore.append(round(dszCode.calculate_similarity_customized(target_struc, product_3_smi),2))
		elif steplength == 4:
			dbidTemp = [reactant_1, product_1, product_2, product_3,product_4]
			compds.append(dbidTemp)
			compdnames.append([getNameFromDBid(i) for i in dbidTemp])

			compdsmis.append([reactant_1_smi, product_1_smi, product_2_smi, product_3_smi,product_4_smi])
			smartsList.append([reactant_1_smarts, product_1_smarts, product_2_smarts, product_3_smarts,product_4_smarts])
			smiScore.append(round(dszCode.calculate_similarity_customized(target_struc, product_4_smi),2))

	single_cluster["compds"] = compds
	single_cluster["compdsmis"] = compdsmis
	single_cluster["compdnames"] = compdnames
	single_cluster["smartsList"] = smartsList
	single_cluster["smiScore"] = smiScore
	single_cluster["bbsRxnID"] = bbsRxnID

	single_cluster = sortedSingleCluster(single_cluster)
	result_total_forone["single_cluster"] = single_cluster

	return result_total_forone
def extendbbsRxnIDs(bbsRxnIDs):
	bbsretroids= list(set([bbs_rxn_merge.objects.filter(bbsRxnID = i).first().bbsretroid for i in bbsRxnIDs]))
	objtemps = [bbs_rxn_merge.objects.filter(bbsretroid=i) for i in bbsretroids]
	objs = []
	bbsRxnIDsNew = []
	for each_one in objtemps:
		objs.extend(each_one)
	return objs
def getSimilaritySolution(request):
	similarity_compd = request.GET.get("similarity_compd").replace("\\",'').replace("/",'')
	target_struc = request.GET.get("target_struc")

	bbsrxnid_bbsretroid = [[i.bbsRxnID,i.bbsretroid] for i in bbs_rxn_merge.objects.filter(endcompdsmi=similarity_compd)]
	bbsretroids = list(set([i[1] for i in bbsrxnid_bbsretroid]))
	bbsRxnIDs = [i[0] for i in bbsrxnid_bbsretroid]
	#bbsRxnIDs = list(set([j.bbsRxnID for i in bbsretroids for j  in bbs_rxn_merge.objects.filter(bbsretroid=i)]))
	# bbsRxnIDs = [i.bbsRxnID for i in bbs_rxn_merge.objects.filter(bbsretroid=i)]
	# bbsRxnIDs = list(set([i[0] for i in bbsrxnid_bbsretroid]))
	#bbsRxnIDs = list(set([i.bbsRxnID for i in bbs_rxn_merge.objects.filter(endcompdsmi=similarity_compd)]))
	#bbsRxnIDs = [18039]
	#newobjs = extendbbsRxnIDs(bbsRxnIDs)
	#bbsRxnIDs.append(202522)
	bbsRxnIDsObjs = [pathsearch_combo_steps_smarts.objects.filter(bbsRxnID = i).first() for i in bbsRxnIDs]
	dict_bbsRxnIDs_count = {i.bbsRxnID:i.combo_simple_smarts_count for i in bbsRxnIDsObjs if i.combo_simple_smarts_count != 0}
	top20 = sorted(dict_bbsRxnIDs_count.items(),key=lambda x:x[1],reverse=True)
	top20 = top20[:20]


	result_total = []

	for each_one in top20:
		bbsRxnID,count = each_one
		smarts = bbs_rxn_merge.objects.get(bbsRxnID = bbsRxnID).smarts
		if smarts == None or smarts == '' or smarts == '>>':
			continue
		precursor = getPredictPrecursorFromTarget(smarts,target_struc)
		result_total.extend(getKnownRouteFrombbsRxnIDs([bbsRxnID],precursor,target_struc))
	
	result_total_new = []
	for i in result_total:
		i["Score"] = round(math.log(i["combo_simple_smarts_count"])+(10-sum(i["SaScore"])/i["steplength"])+(4-i["steplength"]),2) 
		result_total_new.append(i)
	result_total = result_total_new
	result_total1 = sorted([i for i in result_total if i["steplength"] == 1],key=lambda x:x["Score"],reverse=True)
	
	
	result_total2 = sorted([i for i in result_total if i["steplength"] == 2],key=lambda x:x["Score"],reverse=True)
	result_total3 = sorted([i for i in result_total if i["steplength"] == 3],key=lambda x:x["Score"],reverse=True)
	
	result_total4 = sorted([i for i in result_total if i["steplength"] == 4],key=lambda x:x["Score"],reverse=True)
	
	for each_one in result_total1:
		each_one["single_cluster"] = sortedSingleCluster(each_one["single_cluster"])
	for each_one in result_total2:
		each_one["single_cluster"] = sortedSingleCluster(each_one["single_cluster"])
	for each_one in result_total3:
		each_one["single_cluster"] = sortedSingleCluster(each_one["single_cluster"])
	for each_one in result_total4:
		each_one["single_cluster"] = sortedSingleCluster(each_one["single_cluster"])

	result_total = []

	result_total1_new = merge_result_total_n(result_total1)
	result_total1 = [getExtendResult_total(i,target_struc) for i in result_total1_new]

	result_total.extend(result_total1)

	result_total2_new = merge_result_total_n(result_total2)
	result_total2 = [getExtendResult_total(i,target_struc) for i in result_total2_new]

	result_total.extend(result_total2)

	result_total3_new = merge_result_total_n(result_total3)
	result_total3 = [getExtendResult_total(i, target_struc) for i in result_total3_new]


	result_total.extend(result_total3)

	result_total4_new = merge_result_total_n(result_total4)
	result_total4 = [getExtendResult_total(i, target_struc) for i in result_total4_new]

	result_total.extend(result_total4)

	result_total = sorted(result_total,key=lambda x:(x["steplength"],-x['Score']))
	return render(request,"rxncluster_similarity_pathdetail.html",{"result_total":result_total})

def merge_result_total_n(result_totaln):
	result_total_new = []
	product_pre_lists = []
	for each_one in result_totaln:
		product_pre_list = each_one["product_pre_list"]
		single_cluster = each_one["single_cluster"]
		if product_pre_list not in product_pre_lists:
			product_pre_lists.append(product_pre_list)
			result_total_new.append(each_one)
		else:
			combo_simple_smarts_count = each_one["combo_simple_smarts_count"]
			final_each_one = [i for i in result_total_new if i["product_pre_list"]==product_pre_list][0]
			final_each_one['combo_simple_smarts_count'] = final_each_one["combo_simple_smarts_count"]+combo_simple_smarts_count
			single_cluster_final = final_each_one["single_cluster"]
			single_cluster_final["smiScore"].extend(single_cluster["smiScore"])
			single_cluster_final["compdsmis"].extend(single_cluster["compdsmis"])
			single_cluster_final["smartsList"].extend(single_cluster["smartsList"])
			single_cluster_final["compds"].extend(single_cluster["compds"])
			single_cluster_final["compdnames"].extend(single_cluster["compdnames"])
			single_cluster_final["bbsRxnID"].extend(single_cluster["bbsRxnID"])
			final_each_one["single_cluster"] = single_cluster_final
			final_each_one["detail"].extend(each_one["detail"])
	return 	result_total_new

def API_result_request(request,precursor,target):
	new_dict = getResultPrecursorTarget(precursor,target)
	return HttpResponse(json.dumps(new_dict),content_type='application/json')

def API_rule_extract(request,reactants,products):
	smirks = reactants +">>"+products
	result = getRulesFromSmirks(smirks)
	return HttpResponse(json.dumps({"AAM":result[0],"rule":result[2]}),content_type='application/json')




