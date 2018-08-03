from dendropy import Tree
from copy import copy

class TrplDistr(object):
    def __init__(self, ddpTree = None, tree_file = None, schema = "newick"):
        if tree_file:
            self.ddpTree = Tree.get_from_path(tree_file,schema,preserve_underscores=True)
        else:
            self.ddpTree = ddpTree

    def get_all_distributions(self,k):
        for i in range(1,k+1):
            self.compute_distribution(i)

    def add_distributions(self,d1,d2):
        d = copy(d1)
        for (s,c) in d2:
            self.update_distribution(d,s,c)
        return d

    def update_distribution(self,d,score,count):
        for i,(s,c) in enumerate(d):
            if s == score:
                d[i] = (s,c+count)
                return
        d.append((score,count))       

    def compute_distribution(self,k):
# assume we have all the distributions of 1,2,...,k-1
# assume ddpTree has no polytomy
        if k == 1:
            for node in self.ddpTree.postorder_node_iter():
                if node.is_leaf():
                    node.distributions = {1:[(0,1)]}
                else:
                    node.distributions = {1:[(0,sum([c.distributions[1][0][1] for c in node.child_node_iter()]))]} 
        else:
            for node in self.ddpTree.postorder_node_iter():
                if not node.is_leaf():
                    c1,c2 = node.child_nodes()
                    n1 = c1.distributions[1][0][1]
                    n2 = c2.distributions[1][0][1]
                     
                    if n1+n2 < k:
                        continue           
                    
                    if n1 > n2:
                        # swap
                        n = n1
                        n1 = n2
                        n2 = n
                        c = c1
                        c1 = c2
                        c2 = c
                    
                    node.distributions[k] = []

                    if n1 >= k:
                        node.distributions[k] = self.add_distributions(node.distributions[k],c1.distributions[k])
                    if n2 >= k:
                        node.distributions[k] = self.add_distributions(node.distributions[k],c2.distributions[k])

                    m = min(n1,k-1)
                    n3 = self.ddpTree.seed_node.distributions[1][0][1] - n1 - n2
                    for i in range(m,0,-1):
                        j = k-i
                        if i > n1 or j > n2:
                            continue
                        for score_1,count_1 in c1.distributions[i]:
                            for score_2,count_2 in c2.distributions[j]:
                                score = score_1 + score_2 + i*j*n3
                                count = count_1*count_2
                                self.update_distribution(node.distributions[k],score,count)
                                #node.distributions[k].append( (score,count) )

    def report_pdf(self,k):
        f = sorted(self.ddpTree.seed_node.distributions[k])
        s = sum([y for (x,y) in f])
        return [(x,float(y)/s) for (x,y) in f]

    def report_cdf(self,k):
        p = self.report_pdf(k)
        c = [0]*len(p)
        v = 0
        for i,(s,f) in enumerate(p):
            c[i] = (s,v) 
            v += f
        return c    

    
    def report_quantile(self,trplScore,k):
        def __binary_search__(l=0,r=None):
            
            if r is None:
                r = len(c)-1
            med = int((r+l)/2)
            s,q = c[med]

            if s >= trplScore:
                if med == l:
                    return med
                return __binary_search__(l,med-1)
            else:
                if med == r:
                    return r+1
                return __binary_search__(med+1,r) 
        
        c = self.report_cdf(k)
        idx = __binary_search__()
        return c[idx][1]

    def report_pval(self,trplScore,k):
        return 1-self.report_quantile(trplScore,k)
        

def main():
    from sys import argv

    D = TrplDistr(tree_file = argv[1])

    k=int(argv[2])

    D.get_all_distributions(k)

    print(D.report_pdf(k))
    print(D.report_cdf(k))
    print(D.report_quantile(int(argv[3]),k))
    print(D.report_pval(int(argv[3]),k))

if __name__=="__main__":
    main()    
