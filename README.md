# Lagrange_polinomial_interpolation
## Predicts the approximate location of points without an explicit function based on a dataset of known points
n == points data steps 

another_n == Equal x steps to determine the approximate curve by means of the Lagrange polynomial

a == start, b == end of curve

pointsData == vector of {{x_i, y_i},...}.

To set the array manually can be done in this way:
> a=1; 
> b=2; 
> pointsData = {{1.1, 1.2}, {1.2, 1.5}, {1.6, 0.9}, {1.8, -0.2}};

PointsWeKnow_FromSteps() set the pointsData;

offset() - offset from the 'a' point + i * step(equals to (b-a) * 1/n)

e.g.

>![изображение](https://user-images.githubusercontent.com/69731829/133691626-cacd7e61-2b8b-4f30-ac17-55eeeb62d689.png)
>![изображение](https://user-images.githubusercontent.com/69731829/133691394-efd9025e-b3fd-4b9d-bf8a-ffbbc7fd16e8.png)
>![изображение](https://user-images.githubusercontent.com/69731829/133692927-0f1c864f-6606-4b12-b045-f5e07ba46395.png)

