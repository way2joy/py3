# -*- coding: utf-8 -*-
"""
Created on Mon Jul  2 14:15:00 2018

@author: unhochang
"""

#%%
import numpy as np

def rpkm(counts, lengths):
    """RPKM을 계산한다.
    RPKM = (10^9 * C) / (N * L)

    변수 :
    C = 유전자에 매핑된 판독 수
    N = 실험에서 매핑된 총 판독 수
    L = 유전자 염기쌍 엑손(Exon) 길이

    매개변수
    ----------
    counts: array, shape (N_genes, N_samples)
        RNA 염기서열분석 개수 (열 : 개별 샘플, 행 : 유전자)
    lengths: array, shape (N_genes,)
        유전자 행과 같은 순서로 된 염기쌍 유전자 길이

    반환값
    -------
    normed : array, shape (N_genes, N_samples)
        정규화된 RPKM 개수 행렬
    """
    N = np.sum(counts, axis=0)  # 각 열의 합계 (샘플 당 총 판독수)
    L = lengths
    C = counts

    normed = 1e9 * C / (N[np.newaxis, :] * L[:, np.newaxis])

    return(normed)
    
#%%
gene0 = [100, 200]
gene1 = [50, 0]
gene2 = [350, 100]

expression_data = [gene0, gene1, gene2]
expression_data[2][0]

#%%
type(gene0)
type(expression_data)

#%%
import numpy as np

array1d = np.array([1, 2, 3, 4])
print(array1d)
print(type(array1d))
print(array1d.shape)
array1d.ndim

#%%
array2d = np.array(expression_data)
print(array2d)
print(array2d.shape)
print(type(array2d))
print(array2d.ndim) 

#%%
# 정수 범위(0~999,999)의 ndarray 생성 
array = np.arange(1e6)
array.shape
# 파이썬 리스트로 변환
list_array = array.tolist()

%timeit -n10 y = [val * 5 for val in list_array]

%timeit -n10 x = array * 5

#%%
# ndarray x 생성
x = np.array([1, 2, 3], np.int32)
print(x)
print(x.shape)
#%% slicing
y = x[:2]
print(y)

#%%
y[0] = 5
y

#%% 참조와 복사의 차이를 이해해 둘 것.
y = np.copy(x[:2])
y     # y를 바꾸면 x도 같이 바뀐다.

y[0] = 4
print(y)
print(x)

#%%
x = np.array([1,2,3,4])
print(x * 2)
y = np.array([0,1,2,1])
print(x + y)
x
x = np.reshape(x, (len(x), 1))
x
y
y = np.reshape(y, (1, len(y)))
y
x.shape
y.shape

outer = x * y
outer
outer.shape

#%%
import pandas as pd
filename = "data/counts.txt"
with open(filename, 'rt') as f:
    data_table = pd.read_csv(f, index_col=0)

#%%
dir(data_table)
type(data_table)
data_table.iloc[:2,:2]
samples = list(data_table.columns)
len(samples)
samples[:6]

#%%
with open("data/genes.csv", "rt") as f:
    gene_info = pd.read_csv(f, index_col=0)

#%%
type(gene_info)
gene_info.iloc[:5,:]

data_table.shape[0]
gene_info.shape[0]

matched_index = pd.Index.intersection(data_table.index, gene_info.index)
matched_index.shape
matched_index[:10]
data_table.index[:10]

#%%
counts = np.asarray(data_table.loc[matched_index], dtype=int)
type(counts)
counts.shape
counts[:5][:5]

gene_names = np.array(matched_index)
len(gene_names)

# 포맷에 맞춘 프린트 방법을 찾음.
# inline 명령어 실행에 컬리 브레이스 활용.
print(f'{counts.shape[1]}개의 개체에 {counts.shape[0]}개의 유전자 측정됨')

#%%
gene_lengths = np.asarray(gene_info.loc[matched_index]['GeneLength'], dtype=int)
counts.shape
gene_lengths.shape

#%%
import matplotlib.pyplot as plt
plt.style.use('ggplot')

import sys
sys.stdin.encoding
sys.stdout.encoding

#%%
total_counts = np.sum(counts, axis=0)
len(total_counts)

from scipy import stats
density = stats.kde.gaussian_kde(total_counts)

x = np.arange(min(total_counts), max(total_counts), 10000)

#%%
fig, ax = plt.subplots()
ax.plot(x, density(x))
ax.set_xlabel("개체당 유전자 발현 횟수")
ax.set_ylabel("밀도")
## plt.show()

#%%
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm

font_list = fm.findSystemFonts(fontpaths=None, fontext='ttf')
len(font_list)

print(fm.fontManager.ttflist[1:10].index)
dir(fm.fontManager.ttflist)

f = [f.name for f in fm.fontManager.ttflist]
f[:10]
[(f.name, f.fname) for f in fm.fontManager.ttflist if 'Nanum' in f.name]


#%%
# 그래프에서 마이너스 폰트 깨지는 문제에 대한 대처
mpl.rcParams['axes.unicode_minus'] = False
plt.rcParams["font.family"] = 'NanumGothicCoding'
plt.rcParams["font.size"] = 15
plt.rcParams["figure.figsize"] = (10,5)

print('# 설정되어있는 폰트 사이즈')
print (plt.rcParams['font.size'] ) 
print('# 설정되어있는 폰트 글꼴')
print (plt.rcParams['font.family'] )

fig, ax = plt.subplots()
ax.plot(x, density(x))
ax.set_xlabel("개체당 유전자 발현 횟수")
ax.set_ylabel("밀도")
fig
ax

#%%
np.random.seed(seed=7)
sample_index = np.random.choice(range(counts.shape[1]), size=70,
                                replace=False)
sample_index
counts_subset = counts[:,sample_index]

#%%
def reduce_xaxis_label(ax, factor):
    plt.setp(ax.get_ticklabels(), visible=False)
    for label in ax.get_ticklabels()[factor-1::factor]:
        label.set_visible(True)

#%%
plt.boxplot(np.log(counts_subset + 1))

#%%
 
