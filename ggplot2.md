# Data visualization with ggplot2
`ggplot2` is a powerful and flexible R package for creating stunning visualizations. It implements the grammar of graphics, an easy to understand and implement concept in data visualization.

This tutorial will use ggplot2 in R-Studio.


#### references
https://ggplot2.tidyverse.org/reference/

https://ggplot2.tidyverse.org/articles/ggplot2-specs.html

https://r-graphics.org/



## Installation
Install `ggplot2` from CRAN as follows:

```r
install.packages("ggplot2")
```

alternatively, you can use `ggplot2` from the `tidyverse` package.

```r
install.packages("tidyverse")
```
-----
## Basic Usage

To use `ggplot2`, you first need to load it into your R environment. You can do this using the `library()` function:

```r
library(ggplot2)
```

The main function in `ggplot2` is `ggplot()`, which serves as a base layer for all other additions to the plot. 

```r
ggplot(data = <DATA>, mapping = aes(<MAPPINGS>)) +  <GEOM_FUNCTION>()
```

Data: The raw data that you want to plot.

Aesthetics `aes()`: Aesthetics of the geometric and statistical objects, such as position, color, size, shape, and transparency.

Geometries `geom_`: The geometric shapes that will represent the data.

Here's an example:

```r
ggplot(data = mtcars, mapping = aes(x = mpg, y = disp)) + geom_point()
```

In this example, `mtcars` is the dataset, `mpg` and `disp` are variables in the dataset, and `geom_point()` adds a layer of points to the plot.

-----
## Overview of the mtcars dataset
Explore `midwest` and `mtcars` datasets. ggplot uses data that are in long format.

```r
midwest
```


```r 
head(mtcars)
```

The data was extracted from the 1974 Motor Trend US magazine, and comprises fuel consumption and 10 aspects of automobile design and performance for 32 automobiles (1973–74 models).

mpg - Miles per Gallon

cyl - Number of cylinders

disp - displacement, in cubic inches

hp - horsepower

drat - Rear axle ratio

wt - weight

qsec - 1/4 mile time; a measure of acceleration

vs - engine shape

am - transmission; auto or manual

gear - Number of gears

carb - Number of carburetors.




#### Scatter plot with R base graphics
```r
plot(mtcars$wt, mtcars$mpg)
```


#### Basic scatter plot using ggplot
```r
ggplot(data = mtcars, mapping = aes(x = wt, y = mpg)) +
  geom_point()
```
`geom_point()` tells ggplot to add a layer of points to the plot.

Anything you put in the `ggplot()` function can be seen by any geom layers that you add (i.e., these are universal plot settings). This includes the x- and y-axis you set up in `aes()`.

You can also specify aesthetics for a given geom independently of the aesthetics defined globally in the `ggplot()` function.

The `+` sign used to add layers must be placed at the end of each line containing a layer. If, instead, the `+` sign is added in the line before the other layer, ggplot2 will not add the new layer and will return an error message.


#### Adding layers
```r
ggplot(data = mtcars, mapping = aes(x = wt, y = mpg)) +
  geom_point() +
  geom_smooth()
```

#### Customizing aesthetics
```r
ggplot(data = mtcars, mapping = aes(x = wt, y = mpg, color = cyl)) +
  geom_point()
```

```r
ggplot(mtcars, aes(x = wt, y = mpg, color = vs)) +
  geom_point()
```



## Customization

`ggplot2` allows for extensive customization. For example, you can add a title to your plot and label your axes as follows:

```r
ggplot(data = mtcars, aes(x = mpg, y = disp)) +
  geom_point() +
  labs(title = "Miles Per Gallon vs. Displacement",
       x = "Miles Per Gallon",
       y = "Displacement")
```

#### Faceting
```r
ggplot(mtcars, aes(x = mpg, y = disp)) +
  geom_point() +
  facet_wrap(~cyl)
```

#### Customizing scales
```r
ggplot(mtcars, aes(x = mpg, y = disp)) +
  geom_point() +
  scale_x_continuous(limits = c(10, 35))
```

#### Adding statistical transformations
```r
ggplot(mtcars, aes(x = mpg)) +
  geom_histogram(binwidth = 2)
```

#### Themes
```r
ggplot(mtcars, aes(x = mpg, y = disp)) +
  geom_point() +
  theme_minimal()
```


#### Adding Text Labels to Points
```r
ggplot(mtcars, aes(x = mpg, y = disp)) +
  geom_point() +
  geom_text(aes(label = rownames(mtcars)), vjust = 1, hjust = 1)
```

#### Creating a Bar Plot
```r
ggplot(mtcars, aes(x = factor(cyl))) +
  geom_bar() +
  labs(x = "Number of Cylinders", y = "Count")
```

#### Creating a Box Plot with Notches

```r
ggplot(mtcars, aes(x = factor(cyl), y = mpg)) +
  geom_boxplot(notch = TRUE) +
  labs(x = "Number of Cylinders", y = "Miles Per Gallon")
```

#### Creating a Violin Plot

```r
ggplot(mtcars, aes(x = factor(cyl), y = mpg)) +
  geom_violin() +
  labs(x = "Number of Cylinders", y = "Miles Per Gallon")
```

```r
ggplot(mtcars, aes(x = factor(cyl), y = mpg)) +
  geom_violin() +
  geom_jitter(height = 0, width = 0.1) +
  labs(x = "Number of Cylinders", y = "Miles Per Gallon")
```
