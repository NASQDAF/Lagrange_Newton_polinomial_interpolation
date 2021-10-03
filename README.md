# Lagrange&Newton_polinomial_interpolation
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

Offset == a class that stores functions that distribute the segments into which the function is divided. Offset is the part of the dataset that is responsible for the X offset.

e.g.

>![изображение](https://user-images.githubusercontent.com/69731829/133691626-cacd7e61-2b8b-4f30-ac17-55eeeb62d689.png)
>![изображение](https://user-images.githubusercontent.com/69731829/133691394-efd9025e-b3fd-4b9d-bf8a-ffbbc7fd16e8.png)
>Here we have 25 points with some [x]-step at range from [a] to [b]. Initially we calculate by the same step algorithm for 10 points of equal dimension on a dataset of 25 points. We compare it to the required function and obtain the residual value. This is what happens with the Lagrange and Newton Forward and Backward algorithms. Then we compute algorithms by the second step algorithm with the recalculation of points by the main function initially there are 10 points.
>
>![изображение](https://user-images.githubusercontent.com/69731829/135759456-bff813ba-5239-4f03-ae1b-557da637f284.png)



