#coding:utf-8
"""
multiple knock down module
"""
import re
from copy import deepcopy
from cobra.manipulation import undelete_model_genes 
# fluxを調節する関数
 
def change_flux(cobra_model, gene_list=None,
                cumulative_deletions=False,
                disable_orphans=False, wt_flux=None,
                fold_change=None):
    """
    #delete_model_genesを改変したもの。
    model内のfluxを確認して、割合で変化を持たせる。
    cobra_model: モデルを入れる(ex:test.create_test_model())
    gene_list: gene名をリストで渡す(ex:["b1234"])。引数が無ければすべてのgeneでsingleの行われる。doubleにしたい場合は[[リスト],[リスト]]。
    cumulative_deletions: False or True.  If True then any previous deletions will be maintained in the model.
    wt_flux: defoltはNone。wtのfluxと比較したい場合に入れる。(wt_flux.solution.x_dict)
    fold_change: fluxを変更する割合。{gene: value}

    12/31
    fold_changeをdictに。
    薬剤のdoseの違いを再現するため(未完)
    fold_change -> fold_dict

    """
    print "change flux", gene_list
    if not hasattr(cobra_model, '_trimmed') or not cumulative_deletions:
        cobra_model._trimmed = False
        cobra_model._trimmed_genes = []
        cobra_model._trimmed_reactions = {} #Store the old bounds in here.
    spontaneous_re = re.compile('(^|(?<=( |\()))s0001(?=( |\)|$))')
    #Allow a single gene to be fed in as a string instead of a list.
    if not hasattr(gene_list, '__iter__') or \
           hasattr(gene_list, 'id'):  #cobra.Gene has __iter__
        gene_list = [gene_list]
 
    if not hasattr(gene_list[0], 'id'):
        if gene_list[0] in cobra_model.genes:
            tmp_gene_dict = dict([(x.id, x) for x in cobra_model.genes])
        else:
            #assume we're dealing with names if no match to an id
            tmp_gene_dict = dict([(x.name, x) for x in cobra_model.genes])
        gene_list = [tmp_gene_dict[x] for x in gene_list]
 
    # Make the genes non-functional
    [setattr(x, 'functional', False) for x in gene_list]
    the_reactions = set()
    [the_reactions.update(x._reaction) for x in gene_list]
    
    # fold change
    from collections import defaultdict
    fold_dict = defaultdict(lambda: 1.0)
    for k in gene_list:
        for v in k._reaction:
            if fold_dict[v.id] > fold_change[k.id]:
                fold_dict[v.id] = fold_change[k.id]
    del defaultdict
    
    for the_reaction in the_reactions:
        old_lower_bound = the_reaction.lower_bound
        old_upper_bound = the_reaction.upper_bound
        the_gene_reaction_relation = deepcopy(the_reaction.gene_reaction_rule)
        print the_reaction.id
        print the_gene_reaction_relation
        for the_gene in the_reaction._genes:
            the_gene_re = re.compile('(^|(?<=( |\()))%s(?=( |\)|$))'%re.escape(the_gene.id))
            if the_gene in gene_list:
                the_gene_reaction_relation = the_gene_re.sub('False', the_gene_reaction_relation)
            else:
                the_gene_reaction_relation = the_gene_re.sub('True', the_gene_reaction_relation)
        the_gene_reaction_relation = spontaneous_re.sub('True', the_gene_reaction_relation)
        if not eval(the_gene_reaction_relation):
            cobra_model._trimmed_reactions[the_reaction] = (old_lower_bound,
                                                            old_upper_bound)
            # wt_reactionのパターンによって、boundaryを変更
            wt_reaction_flux = wt_flux[the_reaction.id]
            if wt_reaction_flux > 0:
                the_reaction.upper_bound = wt_reaction_flux * fold_dict[the_reaction.id]
            elif wt_reaction_flux < 0:
                the_reaction.lower_bound = wt_reaction_flux * fold_dict[the_reaction.id]
            else:
                pass
            cobra_model._trimmed = True
 
        cobra_model._trimmed_genes =  list(set(cobra_model._trimmed_genes + gene_list))
 
 
 
# fluxの値をだす関数
def multi_knockdown(cobra_model, elements=None,
               method='fba', the_problem='reuse', solver='glpk',
               error_reporting=None, wt_flux=None):# 最後にwt_fluxを追加
    """
    Performs optimization simulations to realize the objective defined
    from cobra_model.reactions[:].objective_coefficients after deleting each gene in
    gene_list from the model.
    
    cobra_model: a cobra.Model object
 
    elements: None or dict. {gene: fold_change} 
 
    method: 'fba' or 'moma'
 
    the_problem: Is None or 'reuse'.
 
    solver: 'glpk', 'gurobi', or 'cplex'. デフォルトの'cglpk'は新バージョン用。
 
    Returns a list of dictionaries: growth_rate_dict, solution_status_dict,
    problem_dict where the key corresponds to each reaction in reaction_list.

    wt_flux: 野生型のfluxを入れる。(ex: solution.x_dict)

    12/28
    element_list to elements. 
    ディクショナリでもらうようにして、なかでelement_list と fc_dictを作成
 
    """
    # 引数の整形
    from copy import copy
    element_list = copy(elements.keys()) # geneのlist
    fc_dict = elements.copy() # fold_changeの値のdict
    del copy

    wt_model = cobra_model.copy() #Original wild-type (wt) model.
    wt_model.id = 'Wild-Type'
    # MOMA constructs combined quadratic models thus we cannot reuse a model
    # generated by the cobra_model.optimize call
    if method.lower() == 'moma':
        the_problem = 'return'
        mutant_model = wt_model.copy() # Need a second model for moma
    else:
        mutant_model = cobra_model
    discard_problems = False
    if the_problem:
        the_problem = 'return'
        discard_problems = True
    the_problem = wt_model.optimize(the_problem=the_problem, solver=solver,
                                       error_reporting=error_reporting)
    wt_f = wt_model.solution.f
    wt_status = wt_model.solution.status
    wt_x = deepcopy(wt_model.solution.x)
    wt_x_dict = deepcopy(wt_model.solution.x_dict)

    if element_list is None:
        element_list = mutant_model.genes

    elif not hasattr(element_list[0], 'id'): # ここに入る
        element_list = map(mutant_model.genes.get_by_id, element_list)
    
    else:
        if mutant_model is not cobra_model:
            element_list = [x.id for x in element_list]
            element_list = map(mutant_model.genes.get_by_id, element_list)
 
 
    wt_problem = the_problem
 
    growth_rate_dict = {}
    solution_status_dict = {}
    problem_dict = {}
    combined_model = None
    flux_dict = {}

    # multi knockdownの場合
    if isinstance(element_list, list):
        print "multiple knockdown simulation >>>"
        
        # fluxを変える関数を実行
        change_flux(mutant_model, element_list, wt_flux=wt_flux, fold_change=fc_dict) # 変異を加える遺伝子
        for the_element in element_list:
            if the_element.id == element_list[0].id:
                mutant_model.id = the_element.id
            else:
                mutant_model.id += "&" + the_element.id

        if mutant_model._trimmed:
            if method.lower() == 'fba':
                the_problem = mutant_model.optimize(the_problem="return",
                                                    solver=solver,
                                                    error_reporting=error_reporting)
                flux_dict[mutant_model.id] = mutant_model.solution.x_dict
                solution_status_dict[mutant_model.id] = mutant_model.solution.status
                growth_rate_dict[mutant_model.id] = mutant_model.solution.f
        
            if discard_problems:
                problem_dict[mutant_model.id] = 'discarded'
            else:
                problem_dict[mutant_model.id] = the_problem
            if not the_problem:
                the_problem = wt_problem
            undelete_model_genes(mutant_model)
       # else just use the wt_f and x
        else:
            if discard_problems:
                problem_dict[mutant_model.id] = 'discarded'
            else:
                problem_dict[mutant_model.id] = wt_problem
            growth_rate_dict[mutant_model.id] = wt_f
            flux_dict[mutant_model.id] = wt_x_dict
            solution_status_dict[mutant_model.id] = wt_status
    

    return(flux_dict, growth_rate_dict, solution_status_dict)

"""
def gene2reaction(genes=None):
    if(genes):
        return None
    else:
        return None
"""

if __name__ == "__main__":
    from cobra.test import create_test_model, ecoli_pickle
    import os, itertools
    
    """
    マルチノックダウン
    """
    
    cobra_model = create_test_model(ecoli_pickle)
    cobra_model.optimize(solver="glpk")
    
    wt_flux = cobra_model.solution.x_dict #WTのflux
    wt_model = cobra_model.copy()
    flux_reaction = []
    target_genes = {"b0149": .3, "b3396": .5}
    
    multi_solution = multi_knockdown(wt_model.copy(), target_genes, wt_flux=wt_flux)
    
    print multi_solution[1]

