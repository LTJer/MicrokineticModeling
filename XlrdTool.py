import xlrd
import numpy as np
import xlsxwriter
import sys
import time



def letter_range(start, stop):
    
    for c in range(ord(start.lower()), ord(stop.lower())+1):
        #yield chr(c)
        yield (c-96)-1
        

def ColNameToNum(name):
    pow = 1
    colNum = 0
    for letter in name[::-1]:
            colNum += (int(letter, 36) -9) * pow
            pow *= 26
    return colNum

def NumToColName(colnum) :
    
    colName = xlsxwriter.utility.xl_col_to_name(colnum-1)
    
    return colName
        
def FetchIndex (excelfile, sheet, val, report = 'both'):
    
    book = xlrd.open_workbook(excelfile, on_demand = True)#.sheet_by_name(sheet)
    workbook = book.sheet_by_name(sheet)
    IndexMatrix_row = np.empty([0,1])
    IndexMatrix_col = np.empty([0,1])
    # start = time.time()
    # print('It took', time.time()-start, 'seconds.')
    for row in range(workbook.nrows): 
        for col in range(workbook.ncols):
            if val in str(workbook.cell_value(row, col)):          
                IndexMatrix_row = np.append(IndexMatrix_row,[[row+1]], axis=0) 
                IndexMatrix_col = np.append(IndexMatrix_col,[[col+1]], axis=0)
    book.release_resources()
    del book
    
    if IndexMatrix_row.size == 0 or IndexMatrix_col.size == 0:
        print(f'No matches found with {val} in {sheet}')
        sys.exit()
        
    if 'co' in report.lower():
        
        if np.all(IndexMatrix_col == IndexMatrix_col[0]):
            # print('Fetch unique column')       
            return int(IndexMatrix_col[0])
        
        else: 
            print(f'Fetch multi-columns at {val}')
            return IndexMatrix_col
        
    elif 'ro' in report.lower():
        
        if np.all(IndexMatrix_row == IndexMatrix_row[0]):
            # print('Fetch unique row')       
            return int(IndexMatrix_row[0])
        
        else: 
            print(f'Fetch multi-rows {val}')
            return IndexMatrix_row
              
    else :
        print('Non-specify format')
        return np.hstack((IndexMatrix_row,IndexMatrix_col))
   
    
    
def IndexFetchValue (excelfile, sheet, rown, coln):
    '''
    support both column number and column name, ex 'AAA'
    * No need for "number-1" to rescale number from 0
    
    '''
    book = xlrd.open_workbook(excelfile)#.sheet_by_name(sheet)
    workbook = book.sheet_by_name(sheet)
    if type(coln) == str : 
       coln =  ColNameToNum(coln)
     
    val = workbook.cell_value(rown-1, coln-1)
    
    book.release_resources()
    del book
    
    return val
    
def WordFetchValue (excelfile, sheet, rowword, colnword):
    
    RowIndex = FetchIndex(excelfile, sheet, rowword, report = 'row')
        
    ColIndex = FetchIndex(excelfile, sheet, colnword, report = 'col')
    
    if type(RowIndex) == int and type(ColIndex) == int :
    
        Val = IndexFetchValue(excelfile, sheet, RowIndex, ColIndex)
        
        # print(rowword, colnword, Val)
        
        return Val
    
    else:
        print('multiple indexes')
       
        
    
        

        
   
    
   
    
#############################debug
# A = FetchIndex('Summary table_v2.xlsx', 'Ea with subsurface', 'S1f', report = 'row')
# print(A) 

# B = FetchValue('Summary table_v2.xlsx', 'Ea with subsurface', 7,'D')
# print(B) 


# rowind = FetchIndex('Summary table_v2.xlsx', 'Ea with subsurface', 'S1f', report = 'row')
# print(rowind)
# colind = FetchIndex('Summary table_v2.xlsx', 'Ea with subsurface', 'Fe', report = 'col')
# print(colind)

# C = FetchValue ('Summary table_v2.xlsx', 'Ea with subsurface', rowind, colind)
# print(C)
# FetchIndex()

# F = WordFetchValue ('Summary table_v2.xlsx', 'Ea with subsurface', 'S1f', 'Fe')
# print(F)


# for row in range(workbook.nrows):
#     for col in range(workbook.ncols):   
#         if 'E1f' in str(workbook.cell_value(row, col)) :
#             print (row, col)


# rowi, rowf = 60, 84
# columni, columnf = 'D', 'U'
# Ea_sheet = np.empty([(rowf-rowi), 0])
 


# for col in letter_range(columni, columnf):
#     Ea_1metal = np.array(workbook.col_values(col, start_rowx=60, end_rowx=84)).reshape((rowf-rowi),1)
#     Ea_sheet = np.append(Ea_sheet, Ea_1metal, axis=1)
#     E1f, E1b, E2f, E2b, E3f, E3b, E4f, E4b, E5f, E7f, E7b, E8f, E8b, E9f, E9b, E10f,\
#         E10b, E11f, E11b, E12f, E13f, E13b, E14f, E14b = [float(i) for i in Ea_1metal]

















