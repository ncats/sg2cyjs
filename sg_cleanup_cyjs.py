# Author: Gergely Zahoranszky-Kohalmi, PhD
#
# Email: gergely.zahoranszky-kohalmi@nih.gov
#
# Organization: National Center for Advancing Translational Sciences (NCATS/NIH)
#
# Script: sg_cleanup_cyjs.py
#
# Aim: Convert the Cytoscape JSON file exported by SmartGraph (https://smartgraph.ncats.io) to proper Cytoscape JSON file.
#
# References: 
#
# [1] Zahoránszky-Kőhalmi, G., Sheils, T. & Oprea, T.I. SmartGraph: a network pharmacology investigation platform.
#     J Cheminform 12, 5 (2020). DOI: https://doi.org/10.1186/s13321-020-0409-9
#
# [2] SmartGraph [https://smartgraph.ncats.io/]
# 
# [3] Cytoscape.js [https://js.cytoscape.org/]
#
# [4] P. Shannondoi et al. Cytoscape: A Software Environment for Integrated Models of Biomolecular Interaction Networks. Genome Res. 2003. 13: 2498-2504
#     DOI: 10.1101/gr.1239303
#
# Ref: https://www.programiz.com/python-programming/json
# Ref: https://networkx.github.io/documentation/networkx-1.10/tutorial/tutorial.html#adding-attributes-to-graphs-nodes-and-edges
# Ref: http://json.parser.online.fr/
# Ref: https://www.geeksforgeeks.org/reading-and-writing-json-to-a-file-in-python/


import json
import sys
import networkx as nx


def stringify_list (l):
	s = ''
	first = True
	
	for i in range(len(l)):
		if first:
			s = l[i]
			first = False
		else:
			s += ';' + l[i]

	return (s)
		


#nodes_target = []
#nodes_compound = []
#nodes_pattern = []

#edges_c2t = []
#edges_p2c = []
#edges_p2t = []
#edges_ppi = []



print ('[SYNTAX] python sg_cleanup_cyjs.py <input_file | exported by SmartGraph to "Cytoscape JSON"> <output_file | use .cyjs  extension!>')

fname = sys.argv[1]
fname_out = sys.argv[2]




with open(fname) as f:
	data = json.load(f)


G = nx.DiGraph()

#print(data['edges'])

edges_json = data['elements']['edges']

nodes_json = data['elements']['nodes']


for n in nodes_json:
	node_record = n['data']
	node_attributes = node_record['node']

	uuid = node_attributes['uuid']

	#print (node_record)

	node_pos = n['position']

	node_type = node_attributes['labels'][0]	

	pos_x = node_pos['x']
	pos_y = node_pos['y']

	G.add_node(uuid)
	G.nodes[uuid]['node_type'] = node_type
	G.nodes[uuid]['pos_x'] = pos_x
	G.nodes[uuid]['y'] = node_pos['y']


	if node_type == 'Target':
		uniprot_id = node_attributes['uniprot_id']
		name = node_attributes['fullname']
		synonyms = node_attributes['synonyms']
		synonyms = stringify_list(synonyms)
		gene_symbol = node_attributes['genes']

		G.nodes[uuid]['uniprot_id'] = uniprot_id
		G.nodes[uuid]['name'] = name
		G.nodes[uuid]['synonyms'] = synonyms
		G.nodes[uuid]['gene_symbol'] = gene_symbol




	elif node_type == 'Compound':
		smiles = node_attributes['hash']
		hash = node_attributes['hash']
		nostereo_hash = node_attributes['nostereo_hash']
		compound_id = node_attributes['compoundId']

		G.nodes[uuid]['smiles'] = smiles
		G.nodes[uuid]['hash'] = hash
		G.nodes[uuid]['nostereo_hash'] = nostereo_hash
		G.nodes[uuid]['copound_id'] = compound_id



		
	elif node_type == 'Pattern':
		smiles = node_attributes['hash']
		hash = node_attributes['hash']
		pattern_id = node_attributes['pattern_id']
		pattern_type = node_attributes['pattern_type']

		G.nodes[uuid]['smiles'] = smiles
		G.nodes[uuid]['hash'] = hash
		G.nodes[uuid]['pattern_id'] = pattern_id
		G.nodes[uuid]['pattern_type'] = pattern_type



	else:
		print ('[ERROR] Invalid node type detected: %s' %(node_type))
		sys.exit(-1)





for e in edges_json:
	edge_record = e['data']
	uuid = edge_record['id']
	
	source_node = edge_record['source']
	target_node = edge_record['target']





	edge_attributes = edge_record['properties']['properties']
	uuid = node_attributes['uuid']
	edge_type = edge_record['properties']['type']
	

	G.add_edge (source_node, target_node)
	G[source_node][target_node]['edge_type'] = edge_type
	
	if edge_type == 'REGULATES':
		edge_info = edge_record['properties']['properties']['edgeInfo']
		edge_info = stringify_list(edge_info)
		source_db = edge_record['properties']['properties']['sourceDB']
		confidence_score = edge_record['properties']['max_confidence_value']
		unique_label = edge_record['properties']['properties']['ppi_uid']
		modulation_type = edge_record['properties']['edgeType']

		G[source_node][target_node]['edge_info'] = edge_info
		G[source_node][target_node]['source_db'] = source_db
		G[source_node][target_node]['confidence_score'] = confidence_score
		G[source_node][target_node]['unique_label'] = unique_label
		G[source_node][target_node]['modulation_type'] = modulation_type
	

	
	elif edge_type == 'TESTED_ON':
		unique_label = edge_record['properties']['properties']['unique_label']
		activity = edge_record['properties']['properties']['activity']
		activity_type = edge_record['properties']['properties']['activity_type']
		modulation_type = edge_record['properties']['edgeType']
	
		G[source_node][target_node]['unique_label'] = unique_label
		G[source_node][target_node]['activity_uM'] = activity
		G[source_node][target_node]['activity_type'] = activity_type
		G[source_node][target_node]['modulation_type'] = modulation_type
	


	

	elif edge_type == 'PATTERN_OF':
		unique_label = edge_record['properties']['properties']['unique_label']
		is_largest = edge_record['properties']['properties']['islargest']
		overlap_ratio = edge_record['properties']['properties']['ratio']

		G[source_node][target_node]['unique_label'] = unique_label
		G[source_node][target_node]['is_largest'] = is_largest
		G[source_node][target_node]['overlap_ratio'] = overlap_ratio
		#G[source_node][target_node]['modulation_type'] = 'NA'
	
	
	elif edge_type == 'POTENT_PATTERN_OF':
		unique_label = edge_record['properties']['properties']['unique_label']
	
		G[source_node][target_node]['unique_label'] = unique_label
		#G[source_node][target_node]['modulation_type'] = 'NA'
	
	else:
		print ('[ERROR] Invalid edge type detected: %s' %(edge_type))
		sys.exit(-1)


	

	
	#G.add_node(uuid)

	#print (node_record)

#print (G)

j = nx.cytoscape_data(G, attrs=None)

with open(fname_out, 'w+') as outfile: 
    json.dump(j, outfile)


print ('[Done.]')

