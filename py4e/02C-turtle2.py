# 거북이 그래픽-도형그리기2
import turtle as t
# 삼각형 그리기
t.color("red")      # <추가> 펜 색상을 빨간색으로 바꿉니다.
t.forward(100)
t.left(120)
t.forward(100)
t.left(120)
t.forward(100)
t.left(120)
# 사각형 그리기
t.color("green")   # <추가> 펜 색상으로 녹색으로 바꿉니다.
t.pensize(3)        # <추가> 펜 두께를 3으로 바꿉니다.
t.forward(100)
t.left(90)
t.forward(100)
t.left(90)
t.forward(100)
t.left(90)
t.forward(100)
t.left(90)
# 원
t.color("blue")     # <추가> 펜 색상을 파란색으로 바꿉니다.
t.pensize(5)        # <추가> 펜 두께를 5로 바꿉니다.
t.circle(50)
