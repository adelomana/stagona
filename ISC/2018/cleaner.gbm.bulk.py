import sys

g=open("GSE57872_GBM_data_matrix.ALO.bulk.txt","w")

with open("GSE57872_GBM_data_matrix.txt", 'r') as f:
    header=f.readline()
    v=header.split('\t')
    v[-1]=v[-1].replace("\n","")

    selected=v[1:]
    
    print(selected,len(selected))

    final=[]
    for element in selected:
        good=False
        for tag in ["MGH26", "MGH28", "MGH29", "MGH30", "MGH31"]:
            if tag in element:
                if "Tumor" in element:
                    final.append(element)
                    break

    print(final,len(final))

    g.write('geneName')
    for element in final:
        finalElement=element.replace("MGH264","MGH26")
        finalElement=finalElement.replace("Tumor","Bulk")
        g.write("\t{}".format(finalElement))
    g.write('\n')
    
    for line in f:
        v=line.split('\t')

        v[-1]=v[-1].replace('\n','')
        
        g.write(v[0])

        for i in range(len(v[1:])):
            
            label=selected[i]
            if label in final:
                g.write("\t{}".format(v[i+1]))        
        g.write("\n")
        
g.close()
