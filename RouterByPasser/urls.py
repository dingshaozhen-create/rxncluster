from django.urls import path,include,re_path
from . import views
import os
from django.views.static import serve
current_path = os.path.dirname(os.path.abspath(__file__))
urlpatterns = [re_path('^$',views.rxnCluster,name="rxnCluster"),
re_path(r'search/$',views.search,name='search'),
re_path(r'smiles_to_mol/$',views.smiles_to_mol,name='smiles_to_mol'),
re_path(r'result/(.*)$',views.result,name='result'),
re_path(r'smi2img/(.*)/(.*)/(.*)$',views.smi2img,name='smi2img'),
re_path(r'retro2img/(.*)/(.*)/(.*)$',views.retro2img,name='retro2img'),
re_path(r'showDetail/$',views.showDetail,name='showDetail'),
re_path(r'marksmiles2img/(.*)/(.*)/(.*)/(.*)$',views.marksmiles2img,name='marksmiles2img'),
re_path(r'getRefRxnFromSmirks/$',views.getRefRxnFromSmirks,name='getRefRxnFromSmirks'),
re_path(r'browse/$',views.browse,name="browse"),
re_path(r'statistics/$',views.statistics,name="statistics"),
re_path(r'faq/$',views.faq,name="faq"),
re_path(r'getSimilaritySolution/$',views.getSimilaritySolution,name='getSimilaritySolution'),
re_path(r'API_result_request/(.*)/(.*)$',views.API_result_request,name='API_result_request'),
re_path(r'API_rule_extract/(.*)/(.*)$',views.API_rule_extract,name='API_rule_extract')

]

urlpatterns.append(re_path(r'static/(?P<path>.*)$',serve,{'document_root':current_path+'/static/'}))
urlpatterns.append(re_path(r'media/(?P<path>.*)$',serve,{'document_root':current_path+'/static/'}))
