# 연속한 숫자 제곱의 합 구하기 알고리즘
# 입력: n
# 출력: 1부터 n까지 연속한 숫자의 제곱을 더한 합

#%%
def sum_real_sq(n):
    return n * (n + 1) * (2 * n + 1) / 6   # // 몫의 정수만 구하는 연산자

#%%
print(sum_real_sq(10))   # 1부터 10까지 제곱의 합(입력:10, 출력:385)
print(sum_real_sq(100))  # 1부터 100까지 제곱의 합(입력:100, 출력:338350)
