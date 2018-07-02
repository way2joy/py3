# 연속한 숫자 제곱의 합을 구하는 알고리즘
# 입력: n
# 출력: 1부터 n까지 연속한 숫자의 제곱을 더한 합

#%%
import os
os.getcwd() # 현재 작업 디렉토리
dir(os)
os.chdir("./algo")

#%%
def gisa_sq(n):
    s = 0
    for i in range(1, n + 1):
        s = s + i * i
    return s

print(gisa_sq(15))  # F9키로 한 줄 시행, Ctrl+엔터로 셀 시행
print(gisa_sq(300))


#%% 정답
def sum_sq(n):
    s = 0
    for i in range(1, n + 1):   
        s = s + i * i
    return s

print(sum_sq(10))   # 1부터 10까지 제곱의 합(입력:10, 출력:385)
print(sum_sq(100))  # 1부터 100까지 제곱의 합(입력:100, 출력:338350)
 
