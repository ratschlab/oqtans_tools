#!/usr/bin/env python2.6
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# Written (W) 2011 Christian Widmer
# Copyright (C) 2011 Max-Planck-Society

"""
Created on 27.05.2011
@author: Christian Widmer
@summary: Deal with hierarchical structure

"""


# yet another graph viz binding


class TreeNode(object):
    """
    Simple graph implementation for hierarchical multitask
    """
    
    
    def __init__(self, name=""):
        """
        define fields
        
        @param name: node name, for leaves this defines the data set identifier
        @type name: str
        """

        self.name = name
        self.pretty_name = name
        self.parent = None
        self.children = []
        self.predictor = None
        self.edge_weight = 0
        self.cost = 1.0


        
    def add_child(self, node, weight=1.0):
        """
        add node as child of current leaf
        
        @param node: child node
        @type node: TreeNode
        """
        
        node.parent = self
        node.edge_weight = weight
        
        self.children.append(node)


    
    def get_data_keys(self):
        """
        fetch dataset names that are connected to leaves below the current node
        this can be used as key to data structures
        
        @return: list of dataset names
        @rtype: list<str>
        """
        
        return [node.name for node in self.get_leaves()]
        
        
        
        
    def get_leaves(self):
        """
        fetch all leaves with breadth first search from node
        
        @return: list of leaves
        @rtype: list<TreeNode>
        """
        
        leaves = []
        
        
        grey_nodes = [self]
        
        while len(grey_nodes)>0:
            
            
            node = grey_nodes.pop(0) #pop first item (go python!)
            
            if len(node.children) == 0:
                leaves.append(node)
            else:
                grey_nodes.extend(node.children)
                
        return leaves
  

    def get_nearest_neighbor(self):
        """
        
        """

        leaves = [leaf for leaf in self.parent.get_leaves() if leaf!=self]

        leftmost = leaves[0]

        return leftmost

   
   
    def get_all_nodes(self):
        """
        fetch all nodes with breadth first search from node
        
        @return: list of nodes
        @rtype: list<TreeNode>
        """
        
        nodes = []
        
        
        grey_nodes = [self]
        
        while len(grey_nodes)>0:
            
            node = grey_nodes.pop(0) #pop first item (go python!)
            nodes.append(node)
            grey_nodes.extend(node.children)
            
                
        return nodes


    def get_node(self, node_name):
        """
        get node from subtree rooted at self by name
        
        @param node_name: name of node to get
        @type node_name: str
        
        @return: node with name node_name
        @rtype: TreeNode
        """
        
        candidates = [node for node in self.get_all_nodes() if node.name==node_name]
        
        assert(len(candidates)==1)
        
        return candidates[0]
        
        

   
    def get_path_root(self):
        """
        fetch all ancesters of current node (excluding self) 
        until root is reached including root  
        
        @return: list of nodes on the path to root
        @rtype: list<TreeNode>
        """
        
                
        nodes_on_path =[]
        
        node = self

        while node != None:                    
            nodes_on_path.append(node)
            node = node.parent

        return nodes_on_path
    
    
    
    def is_root(self):
        """
        returns true if self is the root node  
        
        @return: indicator if self is root
        @rtype: bool
        """
        
        if self.parent == None:
            return True
        
        else:
            return False


    def is_leaf(self):
        """
        returns true if self is a leaf  
        
        @return: indicator if self is root
        @rtype: bool
        """
        
        if len(self.children) == 0:
            return True
        
        else:
            return False

    
    
    def clear_predictors(self):
        """
        removes predictors from all nodes
        """
        
        all_nodes = self.get_all_nodes()
        
        for node in all_nodes:
            node.predictor = None
            
        
        

    def plot(self, file_name="demo", force_num=False, plot_cost=False, plot_B=False):
        """
        visualizes taxonomy with help of the yetanothergraphvizbinding package
        
        a png is created and tried to open with evince (yes, hardcoded for now)   
        
        @return: graph data structure in yapgvb format
        @rtype: yapgvb.Digraph
        """

        import yapgvb
        
        graph = yapgvb.Digraph("my_graph")
    
        #graph.ranksep = 3
        #graph.ratio = "auto"
    
        grey_nodes = [self]
        
        
        counter = 0
        
        name = ""
        if self.name=="" or force_num:
            name = "root" #str(counter) + ": " + self.name
        else:
            name = self.name
            
        
        self.node = graph.add_node(label = name)
        self.node.color = "gray95"
        
        while len(grey_nodes)>0:
           
            node = grey_nodes.pop(0) #pop first item
                
                    
            print node.pretty_name

            #enqueue children
            if node.children != None:
                
                grey_nodes.extend(node.children)
        
                #add edges
                for child_node in node.children:
                    
                    counter += 1

                    child_name = ""
                    if child_node.name=="" or force_num:
                        child_name = str(counter) + ": " + child_node.pretty_name
                    else:
                        child_name = child_node.pretty_name
                    
                    child_node.node = graph.add_node(label = child_name)
                    
                    child_node.node.style = "filled"
                    
                    if child_node.is_leaf():
                        child_node.node.color = "gray80"
                        child_node.node.shape = "doubleoctagon" #"doublecircle"
                        
                    else:
                        child_node.node.color = "gray95"
                    
                    edge = node.node >> child_node.node
                    
                    
                    tmp_label = ""
                    
                    if plot_cost:
                        try:
                            tmp_label += "C=" + str(child_node.cost)
                        except Exception:
                            print "cost attribute not set"
                        
                    if plot_B:
                        tmp_label += "B=" + str(child_node.edge_weight)
                        
                    edge.label = tmp_label


    
        print "Using dot for graph layout..."
        graph.layout(yapgvb.engines.dot)
        #graph.layout(yapgvb.engines.neato)
        #graph.layout(yapgvb.engines.fdp)
        #graph.layout(yapgvb.engines.twopi)
        #graph.layout(yapgvb.engines.circo)
        
        
        demo_formats = [
            yapgvb.formats.png
            #yapgvb.formats.ps
            #yapgvb.formats.xdot
        ]
        
        for myformat in demo_formats:
            filename = file_name + ".%s" % myformat
    
            print "  Rendering .%s ..." % filename
    
            graph.render(filename)            
        
        
        #os.system("evince " + filename + " &")
        
        return graph


def create_flat_hierarchy():
    
    org_list = [
        "Acidobacterium_capsulatum_ATCC_51196_uid59127",
        "Agrobacterium_tumefaciens_C58_uid57865",
        "Bacillus_anthracis_Ames_uid57909",
        "Bacillus_subtilis_168_uid57675",
        "Bacteroides_thetaiotaomicron_VPI_5482_uid62913",
        "Bifidobacterium_longum_NCC2705_uid57939",
        "Clostridium_botulinum_A_ATCC_19397_uid58927",
        "Enterobacter_638_uid58727",
        "Escherichia_coli_BW2952_uid59391",
        "Escherichia_fergusonii_ATCC_35469_uid59375",
        "Flavobacterium_psychrophilum_JIP02_86_uid61627",
        "Fusobacterium_nucleatum_ATCC_25586_uid57885",
        "Klebsiella_pneumoniae_MGH_78578_uid57619",    
        "Mycoplasma_genitalium_G37_uid57707",
        "Mycoplasma_pneumoniae_M129_uid57709",
        "Rhizobium_leguminosarum_bv__viciae_3841_uid57955",
        "Salmonella_enterica_serovar_Heidelberg_SL476_uid58973",
        "Streptococcus_thermophilus_CNRZ1066_uid58221"        
    ]


    root = TreeNode("root")
    
    for org_name in org_list:
        node = TreeNode(org_name)
        root.add_child(node)
        
    #root.plot()

    return root
    
    

def create_tax_four():
    """
    group two species from same genus
    """

    root = TreeNode("root")
    
    escherichia = TreeNode("escherichia")
    root.add_child(escherichia)
    
    bacillus = TreeNode("bacillus")
    root.add_child(bacillus)
    
    e_coli = TreeNode("Escherichia_coli_BW2952_uid59391")
    escherichia.add_child(e_coli)
    e_fergusonii = TreeNode("Escherichia_fergusonii_ATCC_35469_uid59375")
    escherichia.add_child(e_fergusonii)
    
    b_anthracis = TreeNode("Bacillus_anthracis_Ames_uid57909")
    bacillus.add_child(b_anthracis)
    
    b_subtilis = TreeNode("Bacillus_subtilis_168_uid57675")
    bacillus.add_child(b_subtilis)
    
    #root.plot()

    return root

    

def create_tax_four_distant():
    """
    group two species from same genus
    
    #name num_examples
    #Bacillus_anthracis_Ames_uid57909 2004
    #Escherichia_coli_BW2952_uid59391 1366
    #Methanobrevibacter_smithii_ATCC_35061_uid58827 587
    #Sulfolobus_islandicus_M_14_25_uid58849 756
    """

    root = TreeNode("root")
    
    bacteria = TreeNode("bacteria")
    root.add_child(bacteria)
    
    archaea = TreeNode("archaea")
    root.add_child(archaea)
    
    e_coli = TreeNode("Escherichia_coli_BW2952_uid59391")
    e_coli.pretty_name = "e.coli"
    bacteria.add_child(e_coli)
        
    b_anthracis = TreeNode("Bacillus_anthracis_Ames_uid57909")
    b_anthracis.pretty_name = "b.anthracis"
    bacteria.add_child(b_anthracis)
    
    
    m_smithii = TreeNode("Methanobrevibacter_smithii_ATCC_35061_uid58827")
    m_smithii.pretty_name = "m.smithii"
    archaea.add_child(m_smithii)
    
    s_islandicus = TreeNode("Sulfolobus_islandicus_M_14_25_uid58849")
    s_islandicus.pretty_name = "s.islandicus"
    archaea.add_child(s_islandicus)
    
    
    #root.plot()

    return root

def create_tax_four_distant2():
    """
    group two species from same genus
    
    #name num_examples
    #Escherichia_coli_BW2952_uid59391 1366
    #Clostridium_botulinum_A_ATCC_19397_uid58927
    #Methanobrevibacter_smithii_ATCC_35061_uid58827 587
    #Sulfolobus_islandicus_M_14_25_uid58849 756
    """

    root = TreeNode("root")
    
    bacteria = TreeNode("bacteria")
    root.add_child(bacteria)
    
    archaea = TreeNode("archaea")
    root.add_child(archaea)
    
    e_coli = TreeNode("Escherichia_coli_BW2952_uid59391")
    bacteria.add_child(e_coli)
        
    c_botulinum = TreeNode("Clostridium_botulinum_A_ATCC_19397_uid58927")
    bacteria.add_child(c_botulinum)
    
    m_smithii = TreeNode("Methanobrevibacter_smithii_ATCC_35061_uid58827")
    archaea.add_child(m_smithii)
    
    s_islandicus = TreeNode("Sulfolobus_islandicus_M_14_25_uid58849")
    archaea.add_child(s_islandicus)
    
    
    #root.plot()

    return root



def create_tax_four_distant3():
    """
    group two species from same genus
    
    #name num_examples
    #Bacillus_anthracis_Ames_uid57909 2004
    #Clostridium_botulinum_A_ATCC_19397_uid58927
    #Methanobrevibacter_smithii_ATCC_35061_uid58827 587
    #Sulfolobus_islandicus_M_14_25_uid58849 756
    """

    root = TreeNode("root")
    
    bacteria = TreeNode("bacteria")
    root.add_child(bacteria)
    
    archaea = TreeNode("archaea")
    root.add_child(archaea)
    
    b_anthracis = TreeNode("Bacillus_anthracis_Ames_uid57909")
    #b_anthracis.pretty_name = "b.anthracis"
    bacteria.add_child(b_anthracis)
            
    c_botulinum = TreeNode("Clostridium_botulinum_A_ATCC_19397_uid58927")
    bacteria.add_child(c_botulinum)
    
    m_smithii = TreeNode("Methanobrevibacter_smithii_ATCC_35061_uid58827")
    archaea.add_child(m_smithii)
    
    s_islandicus = TreeNode("Sulfolobus_islandicus_M_14_25_uid58849")
    archaea.add_child(s_islandicus)
    
    
    #root.plot()

    return root


def create_tax_four_distant4():
    """
    group two species from same genus
    
    #name num_examples
    #Bifidobacterium_longum_NCC2705_uid57939 671
    #Clostridium_botulinum_A_ATCC_19397_uid58927
    #Methanobrevibacter_smithii_ATCC_35061_uid58827 587
    #Sulfolobus_islandicus_M_14_25_uid58849 756
    """

    root = TreeNode("root")
    
    bacteria = TreeNode("bacteria")
    root.add_child(bacteria)
    
    archaea = TreeNode("archaea")
    root.add_child(archaea)
    
    b_longum = TreeNode("Bifidobacterium_longum_NCC2705_uid57939")
    #b_longum.pretty_name = "B.longum"
    bacteria.add_child(b_longum)
            
    c_botulinum = TreeNode("Clostridium_botulinum_A_ATCC_19397_uid58927")
    bacteria.add_child(c_botulinum)
    
    m_smithii = TreeNode("Methanobrevibacter_smithii_ATCC_35061_uid58827")
    archaea.add_child(m_smithii)
    
    s_islandicus = TreeNode("Sulfolobus_islandicus_M_14_25_uid58849")
    archaea.add_child(s_islandicus)
    
    
    #root.plot()

    return root


def create_tax_four_distant5():
    """
    group two species from same genus
    
    #name num_examples
    #Klebsiella_pneumoniae_MGH_78578_uid57619 1605
    #Streptococcus_thermophilus_CNRZ1066_uid58221 621
    #Methanobrevibacter_smithii_ATCC_35061_uid58827 587
    #Sulfolobus_islandicus_M_14_25_uid58849 756
    """

    root = TreeNode("root")
    
    bacteria = TreeNode("bacteria")
    root.add_child(bacteria)
    
    archaea = TreeNode("archaea")
    root.add_child(archaea)
    
    k_pneumoniae = TreeNode("Klebsiella_pneumoniae_MGH_78578_uid57619")
    #b_longum.pretty_name = "B.longum"
    bacteria.add_child(k_pneumoniae)
    
    s_thermophilus = TreeNode("Streptococcus_thermophilus_CNRZ1066_uid58221")
    bacteria.add_child(s_thermophilus)
   
    m_smithii = TreeNode("Methanobrevibacter_smithii_ATCC_35061_uid58827")
    archaea.add_child(m_smithii)
    
    s_islandicus = TreeNode("Sulfolobus_islandicus_M_14_25_uid58849")
    archaea.add_child(s_islandicus)
    
    
    #root.plot()

    return root


def create_tax_eight2():
    """
    group two species from same genus
    

    """

    root = TreeNode("root")
   
    bacteria = TreeNode("bacteria")
    root.add_child(bacteria)

    archaea = TreeNode("archaea")
    root.add_child(archaea)
 
    proteobacteria = TreeNode("proteobacteria")
    bacteria.add_child(proteobacteria) 
    
    agrobacterium = TreeNode("Agrobacterium_tumefaciens_C58_uid57865")
    agrobacterium.pretty_name = "A.tumefaciens"
    proteobacteria.add_child(agrobacterium)

    helicobakter = TreeNode("Helicobacter_pylori_26695_uid57787")
    helicobakter.pretty_name = "H.pylori"
    proteobacteria.add_child(helicobakter)
    
    bacillus = TreeNode("bacillus")
    bacteria.add_child(bacillus)
    
    b_anthracis = TreeNode("Bacillus_anthracis_Ames_uid57909")
    b_anthracis.pretty_name = "B.anthracis"
    bacillus.add_child(b_anthracis)
    
    b_subtilis = TreeNode("Bacillus_subtilis_168_uid57675")
    b_subtilis.pretty_name = "B.subtilis"
    bacillus.add_child(b_subtilis)

    # firmicutes
    firmicutes = TreeNode("firmicutes")
    bacteria.add_child(firmicutes)
    
    c_botulinum = TreeNode("Clostridium_botulinum_A_ATCC_19397_uid58927")
    c_botulinum.pretty_name = "C.botulinum"
    firmicutes.add_child(c_botulinum)
    
    s_thermophilus = TreeNode("Streptococcus_thermophilus_CNRZ1066_uid58221")
    s_thermophilus.pretty_name = "S.thermophilus"
    firmicutes.add_child(s_thermophilus)
   
     
    m_smithii = TreeNode("Methanobrevibacter_smithii_ATCC_35061_uid58827")
    m_smithii.pretty_name = "M.smithii"
    archaea.add_child(m_smithii)
    
    s_islandicus = TreeNode("Sulfolobus_islandicus_M_14_25_uid58849")
    s_islandicus.pretty_name = "S.islandicus"
    archaea.add_child(s_islandicus)
     
    root.plot()

    return root


def create_tax_eight():
    """
    group two species from same genus
    

    """

    root = TreeNode("root")
   
    bacteria = TreeNode("bacteria")
    root.add_child(bacteria)

    archaea = TreeNode("archaea")
    root.add_child(archaea)
 
    proteobacteria = TreeNode("proteobacteria")
    bacteria.add_child(proteobacteria) 
    
    escherichia = TreeNode("escherichia")
    proteobacteria.add_child(escherichia)

    agrobacterium = TreeNode("Agrobacterium_tumefaciens_C58_uid57865")
    proteobacteria.add_child(agrobacterium)

    helicobakter = TreeNode("Helicobacter_pylori_26695_uid57787")
    proteobacteria.add_child(helicobakter)

    e_coli = TreeNode("Escherichia_coli_BW2952_uid59391")
    escherichia.add_child(e_coli)
    
    e_fergusonii = TreeNode("Escherichia_fergusonii_ATCC_35469_uid59375")
    escherichia.add_child(e_fergusonii)
    
    bacillus = TreeNode("bacillus")
    bacteria.add_child(bacillus)
    
    b_anthracis = TreeNode("Bacillus_anthracis_Ames_uid57909")
    bacillus.add_child(b_anthracis)
    
    b_subtilis = TreeNode("Bacillus_subtilis_168_uid57675")
    bacillus.add_child(b_subtilis)
   
     
    m_smithii = TreeNode("Methanobrevibacter_smithii_ATCC_35061_uid58827")
    archaea.add_child(m_smithii)
    
    s_islandicus = TreeNode("Sulfolobus_islandicus_M_14_25_uid58849")
    archaea.add_child(s_islandicus)
     
    #root.plot()

    return root



def create_tax_many():
    """
    http://www.ncbi.nlm.nih.gov/genomes/lproks.cgi
    
    Helicobacter_pylori_26695_uid57787,                    #Bacteria; Proteobacteria; Epsilonproteobacteria; Campylobacterales; Helicobacteraceae; Helicobacter; Helicobacter pylori; Helicobacter pylori 26695
    
    Agrobacterium_tumefaciens_C58_uid57865,                #Bacteria; Proteobacteria; Alphaproteobacteria; Rhizobiales; Rhizobiaceae; Rhizobium/Agrobacterium group; Agrobacterium; Agrobacterium tumefaciens; Agrobacterium tumefaciens str. C58
    Rhizobium_leguminosarum_bv__viciae_3841_uid57955,      #Bacteria; Proteobacteria; Alphaproteobacteria; Rhizobiales; Rhizobiaceae; Rhizobium/Agrobacterium group; Rhizobium; Rhizobium leguminosarum; Rhizobium leguminosarum bv. viciae 3841
    
    Escherichia_coli_BW2952_uid59391,                      #Bacteria; Proteobacteria; Gammaproteobacteria; Enterobacteriales; Enterobacteriaceae; Escherichia; Escherichia coli; 
    Escherichia_fergusonii_ATCC_35469_uid59375,            #Bacteria; Proteobacteria; Gammaproteobacteria; Enterobacteriales; Enterobacteriaceae; Escherichia; Escherichia fergusonii; Escherichia fergusonii ATCC 35469
    Enterobacter_638_uid58727,                             #Bacteria; Proteobacteria; Gammaproteobacteria; Enterobacteriales; Enterobacteriaceae; Enterobacter; Enterobacter sp. 638
    Klebsiella_pneumoniae_MGH_78578_uid57619,              #Bacteria; Proteobacteria; Gammaproteobacteria; Enterobacteriales; Enterobacteriaceae; Klebsiella; Klebsiella pneumoniae; Klebsiella pneumoniae subsp. pneumoniae MGH 78578
    Salmonella_enterica_serovar_Heidelberg_SL476_uid58973, #Bacteria; Proteobacteria; Gammaproteobacteria; Enterobacteriales; Enterobacteriaceae; Salmonella; Salmonella enterica; Salmonella enterica subsp. enterica serovar Typhi; Salmonella enterica subsp. enterica serovar Typhi str. CT18
    
    
    Mycobacterium_tuberculosis_H37Rv_uid57777,             #Bacteria; Actinobacteria; Actinobacteridae; Actinomycetales; Corynebacterineae; Mycobacteriaceae; Mycobacterium; Mycobacterium tuberculosis complex; Mycobacterium tuberculosis; Mycobacterium tuberculosis H37Rv
    Bifidobacterium_longum_NCC2705_uid57939,               #Bacteria; Actinobacteria; Actinobacteridae; Bifidobacteriales; Bifidobacteriaceae; Bifidobacterium; Bifidobacterium longum; Bifidobacterium longum subsp. infantis ATCC 15697
    
    Bacillus_anthracis_Ames_uid57909,                      #Bacteria; Firmicutes; Bacillales; Bacillaceae; Bacillus; Bacillus cereus group; Bacillus anthracis; Bacillus anthracis str. 'Ames Ancestor'
    Bacillus_subtilis_168_uid57675,                        #Bacteria; Firmicutes; Bacillales; Bacillaceae; Bacillus; Bacillus subtilis; Bacillus subtilis subsp. natto BEST195
    
    Bacteroides_thetaiotaomicron_VPI_5482_uid62913,        #Bacteria; Bacteroidetes; Bacteroidia; Bacteroidales; Bacteroidaceae; Bacteroides; Bacteroides thetaiotaomicron; Bacteroides thetaiotaomicron VPI-5482
    Flavobacterium_psychrophilum_JIP02_86_uid61627,        #Bacteria; Bacteroidetes; Flavobacteria; Flavobacteriales; Flavobacteriaceae; Flavobacterium; Flavobacterium psychrophilum; Flavobacterium psychrophilum JIP02/86
    
    Mycoplasma_genitalium_G37_uid57707,                    #Bacteria; Tenericutes; Mollicutes; Mycoplasmataceae; Mycoplasma; Mycoplasma genitalium; Mycoplasma genitalium G37
    Mycoplasma_pneumoniae_M129_uid57709,                   #Bacteria; Tenericutes; Mollicutes; Mycoplasmataceae; Mycoplasma; Mycoplasma pneumoniae; Mycoplasma pneumoniae M129
    
    Clostridium_botulinum_A_ATCC_19397_uid58927,           #Bacteria; Firmicutes; Clostridia; Clostridiales; Clostridiaceae; Clostridium; Clostridium botulinum; Clostridium botulinum A str. ATCC 3502
    Streptococcus_thermophilus_CNRZ1066_uid58221           #Bacteria; Firmicutes; Lactobacillales; Streptococcaceae; Streptococcus; Streptococcus thermophilus; Streptococcus thermophilus CNRZ1066
            
    Acidobacterium_capsulatum_ATCC_51196_uid59127,         #Bacteria; Acidobacteria; Acidobacteriales; Acidobacteriaceae; Acidobacterium; Acidobacterium capsulatum; Acidobacterium capsulatum ATCC 51196   

    Fusobacterium_nucleatum_ATCC_25586_uid57885,           #Bacteria; Fusobacteria; Fusobacteriales; Fusobacteriaceae; Fusobacterium; Fusobacterium nucleatum; Fusobacterium nucleatum subsp. nucleatum ATCC 25586
    """


    root = TreeNode("root")

    archaea = TreeNode("archaea")
    root.add_child(archaea)
 
    m_smithii = TreeNode("Methanobrevibacter_smithii_ATCC_35061_uid58827")
    archaea.add_child(m_smithii)
    
    s_islandicus = TreeNode("Sulfolobus_islandicus_M_14_25_uid58849")
    archaea.add_child(s_islandicus)
    
    bacteria = TreeNode("bacteria")
    root.add_child(bacteria)
    
    proteobacteria = TreeNode("proteobacteria")
    bacteria.add_child(proteobacteria) 

    # epsilon    
    helicobakter = TreeNode("Helicobacter_pylori_26695_uid57787")
    proteobacteria.add_child(helicobakter)
    
    # alpha
    alphaproteobacteria = TreeNode("alphaproteobacteria")
    proteobacteria.add_child(alphaproteobacteria)
    
    agrobacterium = TreeNode("Agrobacterium_tumefaciens_C58_uid57865")
    alphaproteobacteria.add_child(agrobacterium)

    rhizobium = TreeNode("Rhizobium_leguminosarum_bv__viciae_3841_uid57955")
    alphaproteobacteria.add_child(rhizobium)
    
    # gamma
    gammaproteobacteria = TreeNode("gammaproteobacteria")
    proteobacteria.add_child(gammaproteobacteria)
    
    escherichia = TreeNode("escherichia")
    gammaproteobacteria.add_child(escherichia)

    e_coli = TreeNode("Escherichia_coli_BW2952_uid59391")
    escherichia.add_child(e_coli)
    
    e_fergusonii = TreeNode("Escherichia_fergusonii_ATCC_35469_uid59375")
    escherichia.add_child(e_fergusonii)

    enterobacter = TreeNode("Enterobacter_638_uid58727")
    gammaproteobacteria.add_child(enterobacter)

    klebsiella = TreeNode("Klebsiella_pneumoniae_MGH_78578_uid57619")
    gammaproteobacteria.add_child(klebsiella)

    salmonella = TreeNode("Salmonella_enterica_serovar_Heidelberg_SL476_uid58973")
    gammaproteobacteria.add_child(salmonella)


    # actino
    actinobacteria = TreeNode("actinobacteria")
    bacteria.add_child(actinobacteria)
    
    m_tuberculosis = TreeNode("Mycobacterium_tuberculosis_H37Rv_uid57777")
    actinobacteria.add_child(m_tuberculosis)
    
    b_longum = TreeNode("Bifidobacterium_longum_NCC2705_uid57939")
    actinobacteria.add_child(b_longum)
    
        
    # bacillus
    bacillus = TreeNode("bacillus")
    bacteria.add_child(bacillus)
    
    b_anthracis = TreeNode("Bacillus_anthracis_Ames_uid57909")
    bacillus.add_child(b_anthracis)
    
    b_subtilis = TreeNode("Bacillus_subtilis_168_uid57675")
    bacillus.add_child(b_subtilis)
    
    
    # bacteroidetes
    bacteroidetes = TreeNode("bacteroidetes")
    bacteria.add_child(bacteroidetes)
    
    c_botulinum = TreeNode("Clostridium_botulinum_A_ATCC_19397_uid58927")
    bacteroidetes.add_child(c_botulinum)
    
    f_psychrophilum = TreeNode("Flavobacterium_psychrophilum_JIP02_86_uid61627")
    bacteroidetes.add_child(f_psychrophilum)
    
    
    # firmicutes
    firmicutes = TreeNode("firmicutes")
    bacteria.add_child(firmicutes)
    
    raise Exception("FIX BUG IN TAXONOMY, botulinum in here twice")
    c_botulinum = TreeNode("Clostridium_botulinum_A_ATCC_19397_uid58927")
    firmicutes.add_child(c_botulinum)
    
    s_thermophilus = TreeNode("Streptococcus_thermophilus_CNRZ1066_uid58221")
    firmicutes.add_child(s_thermophilus)

    # remaining stuff 
    a_capsulatum = TreeNode("Acidobacterium_capsulatum_ATCC_51196_uid59127")
    bacteria.add_child(a_capsulatum)

    f_nucleatum = TreeNode("Fusobacterium_nucleatum_ATCC_25586_uid57885")
    bacteria.add_child(f_nucleatum)

    #root.plot()


    return root



def create_tax_all():
    """
    http://www.ncbi.nlm.nih.gov/genomes/lproks.cgi
    
    Helicobacter_pylori_26695_uid57787,                    #Bacteria; Proteobacteria; Epsilonproteobacteria; Campylobacterales; Helicobacteraceae; Helicobacter; Helicobacter pylori; Helicobacter pylori 26695
    
    Agrobacterium_tumefaciens_C58_uid57865,                #Bacteria; Proteobacteria; Alphaproteobacteria; Rhizobiales; Rhizobiaceae; Rhizobium/Agrobacterium group; Agrobacterium; Agrobacterium tumefaciens; Agrobacterium tumefaciens str. C58
    Rhizobium_leguminosarum_bv__viciae_3841_uid57955,      #Bacteria; Proteobacteria; Alphaproteobacteria; Rhizobiales; Rhizobiaceae; Rhizobium/Agrobacterium group; Rhizobium; Rhizobium leguminosarum; Rhizobium leguminosarum bv. viciae 3841
    
    Escherichia_coli_BW2952_uid59391,                      #Bacteria; Proteobacteria; Gammaproteobacteria; Enterobacteriales; Enterobacteriaceae; Escherichia; Escherichia coli; 
    Escherichia_fergusonii_ATCC_35469_uid59375,            #Bacteria; Proteobacteria; Gammaproteobacteria; Enterobacteriales; Enterobacteriaceae; Escherichia; Escherichia fergusonii; Escherichia fergusonii ATCC 35469
    Enterobacter_638_uid58727,                             #Bacteria; Proteobacteria; Gammaproteobacteria; Enterobacteriales; Enterobacteriaceae; Enterobacter; Enterobacter sp. 638
    Klebsiella_pneumoniae_MGH_78578_uid57619,              #Bacteria; Proteobacteria; Gammaproteobacteria; Enterobacteriales; Enterobacteriaceae; Klebsiella; Klebsiella pneumoniae; Klebsiella pneumoniae subsp. pneumoniae MGH 78578
    Salmonella_enterica_serovar_Heidelberg_SL476_uid58973, #Bacteria; Proteobacteria; Gammaproteobacteria; Enterobacteriales; Enterobacteriaceae; Salmonella; Salmonella enterica; Salmonella enterica subsp. enterica serovar Typhi; Salmonella enterica subsp. enterica serovar Typhi str. CT18
    
    
    Mycobacterium_tuberculosis_H37Rv_uid57777,             #Bacteria; Actinobacteria; Actinobacteridae; Actinomycetales; Corynebacterineae; Mycobacteriaceae; Mycobacterium; Mycobacterium tuberculosis complex; Mycobacterium tuberculosis; Mycobacterium tuberculosis H37Rv
    Bifidobacterium_longum_NCC2705_uid57939,               #Bacteria; Actinobacteria; Actinobacteridae; Bifidobacteriales; Bifidobacteriaceae; Bifidobacterium; Bifidobacterium longum; Bifidobacterium longum subsp. infantis ATCC 15697
    
    Bacillus_anthracis_Ames_uid57909,                      #Bacteria; Firmicutes; Bacillales; Bacillaceae; Bacillus; Bacillus cereus group; Bacillus anthracis; Bacillus anthracis str. 'Ames Ancestor'
    Bacillus_subtilis_168_uid57675,                        #Bacteria; Firmicutes; Bacillales; Bacillaceae; Bacillus; Bacillus subtilis; Bacillus subtilis subsp. natto BEST195
    
    Bacteroides_thetaiotaomicron_VPI_5482_uid62913,        #Bacteria; Bacteroidetes; Bacteroidia; Bacteroidales; Bacteroidaceae; Bacteroides; Bacteroides thetaiotaomicron; Bacteroides thetaiotaomicron VPI-5482
    Flavobacterium_psychrophilum_JIP02_86_uid61627,        #Bacteria; Bacteroidetes; Flavobacteria; Flavobacteriales; Flavobacteriaceae; Flavobacterium; Flavobacterium psychrophilum; Flavobacterium psychrophilum JIP02/86
    
    Mycoplasma_genitalium_G37_uid57707,                    #Bacteria; Tenericutes; Mollicutes; Mycoplasmataceae; Mycoplasma; Mycoplasma genitalium; Mycoplasma genitalium G37
    Mycoplasma_pneumoniae_M129_uid57709,                   #Bacteria; Tenericutes; Mollicutes; Mycoplasmataceae; Mycoplasma; Mycoplasma pneumoniae; Mycoplasma pneumoniae M129
    
    Clostridium_botulinum_A_ATCC_19397_uid58927,           #Bacteria; Firmicutes; Clostridia; Clostridiales; Clostridiaceae; Clostridium; Clostridium botulinum; Clostridium botulinum A str. ATCC 3502
    Streptococcus_thermophilus_CNRZ1066_uid58221           #Bacteria; Firmicutes; Lactobacillales; Streptococcaceae; Streptococcus; Streptococcus thermophilus; Streptococcus thermophilus CNRZ1066
            
    Acidobacterium_capsulatum_ATCC_51196_uid59127,         #Bacteria; Acidobacteria; Acidobacteriales; Acidobacteriaceae; Acidobacterium; Acidobacterium capsulatum; Acidobacterium capsulatum ATCC 51196   

    Fusobacterium_nucleatum_ATCC_25586_uid57885,           #Bacteria; Fusobacteria; Fusobacteriales; Fusobacteriaceae; Fusobacterium; Fusobacterium nucleatum; Fusobacterium nucleatum subsp. nucleatum ATCC 25586
    """
    
    
    root = TreeNode("root")
    
    proteobacteria = TreeNode("proteobacteria")
    root.add_child(proteobacteria) 

    # epsilon    
    helicobakter = TreeNode("Helicobacter_pylori_26695_uid57787")
    proteobacteria.add_child(helicobakter)
    
    # alpha
    alphaproteobacteria = TreeNode("alphaproteobacteria")
    proteobacteria.add_child(alphaproteobacteria)
    
    agrobacterium = TreeNode("Agrobacterium_tumefaciens_C58_uid57865")
    alphaproteobacteria.add_child(agrobacterium)

    rhizobium = TreeNode("Rhizobium_leguminosarum_bv__viciae_3841_uid57955")
    alphaproteobacteria.add_child(rhizobium)
    
    # gamma
    gammaproteobacteria = TreeNode("gammaproteobacteria")
    proteobacteria.add_child(gammaproteobacteria)
    
    escherichia = TreeNode("escherichia")
    gammaproteobacteria.add_child(escherichia)

    e_coli = TreeNode("Escherichia_coli_BW2952_uid59391")
    escherichia.add_child(e_coli)
    
    e_fergusonii = TreeNode("Escherichia_fergusonii_ATCC_35469_uid59375")
    escherichia.add_child(e_fergusonii)

    enterobacter = TreeNode("Enterobacter_638_uid58727")
    gammaproteobacteria.add_child(enterobacter)

    klebsiella = TreeNode("Klebsiella_pneumoniae_MGH_78578_uid57619")
    gammaproteobacteria.add_child(klebsiella)

    salmonella = TreeNode("Salmonella_enterica_serovar_Heidelberg_SL476_uid58973")
    gammaproteobacteria.add_child(salmonella)


    # actino
    actinobacteria = TreeNode("actinobacteria")
    root.add_child(actinobacteria)
    
    m_tuberculosis = TreeNode("Mycobacterium_tuberculosis_H37Rv_uid57777")
    actinobacteria.add_child(m_tuberculosis)
    
    b_longum = TreeNode("Bifidobacterium_longum_NCC2705_uid57939")
    actinobacteria.add_child(b_longum)
    
        
    # bacillus
    bacillus = TreeNode("bacillus")
    root.add_child(bacillus)
    
    b_anthracis = TreeNode("Bacillus_anthracis_Ames_uid57909")
    bacillus.add_child(b_anthracis)
    
    b_subtilis = TreeNode("Bacillus_subtilis_168_uid57675")
    bacillus.add_child(b_subtilis)
    
    
    # bacteroidetes
    bacteroidetes = TreeNode("bacteroidetes")
    root.add_child(bacteroidetes)
    
    c_botulinum = TreeNode("Clostridium_botulinum_A_ATCC_19397_uid58927")
    bacteroidetes.add_child(c_botulinum)
    
    f_psychrophilum = TreeNode("Flavobacterium_psychrophilum_JIP02_86_uid61627")
    bacteroidetes.add_child(f_psychrophilum)
    
    
    # mycoplasma
    mycoplasma = TreeNode("mycoplasma")
    root.add_child(mycoplasma)
    
    m_genitalium = TreeNode("Mycoplasma_genitalium_G37_uid57707")
    mycoplasma.add_child(m_genitalium)
    
    m_pneumoniae = TreeNode("Mycoplasma_pneumoniae_M129_uid57709")
    mycoplasma.add_child(m_pneumoniae)


    # firmicutes
    firmicutes = TreeNode("firmicutes")
    root.add_child(firmicutes)
    
    c_botulinum = TreeNode("Clostridium_botulinum_A_ATCC_19397_uid58927")
    firmicutes.add_child(c_botulinum)
    
    s_thermophilus = TreeNode("Streptococcus_thermophilus_CNRZ1066_uid58221")
    firmicutes.add_child(s_thermophilus)

    # remaining stuff 
    a_capsulatum = TreeNode("Acidobacterium_capsulatum_ATCC_51196_uid59127")
    root.add_child(a_capsulatum)

    f_nucleatum = TreeNode("Fusobacterium_nucleatum_ATCC_25586_uid57885")
    root.add_child(f_nucleatum)

    #root.plot()


    return root


def create_example_graph():
    """
    creates and plots small example
    """

    root = TreeNode("root")
    e_coli = TreeNode("e.coli")
    entero_bakter = TreeNode("entero_bakter")

    weight = 1.0
    root.add_child(e_coli, weight)
    root.add_child(entero_bakter, weight)
    
    #root.plot()

    return root



if __name__ == '__main__':
    create_tax_many()
    
    
