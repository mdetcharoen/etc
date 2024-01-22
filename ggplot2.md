# Data visualization with ggplot2
`ggplot2` is a powerful and flexible R package for creating stunning visualizations. It implements the grammar of graphics, an easy to understand and implement concept in data visualization.

This tutorial will use ggplot2 in R-Studio.

## Installation
Install `ggplot2` from CRAN as follows:

```r
install.packages("ggplot2")
```

alternatively, you can use `ggplot2` from the `tidyverse` package.

```r
install.packages("tidyverse")
```

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

Aesthetics aes(): Aesthetics of the geometric and statistical objects, such as position, color, size, shape, and transparency.

Geometries geom_: The geometric shapes that will represent the data.

Here's an example:

```r
ggplot(data = mtcars, aes(x = mpg, y = disp)) + geom_point()
```

In this example, `mtcars` is the dataset, `mpg` and `disp` are variables in the dataset, and `geom_point()` adds a layer of points to the plot.

### Overview of the mtcars dataset
```r 
head(mtcars)
```

#### Basic scatter plot
```r
ggplot(mtcars, aes(x = mpg, y = disp)) +
  geom_point()
```

#### Adding layers
```r
ggplot(mtcars, aes(x = mpg, y = disp)) +
  geom_point() +
  geom_smooth()
```

#### Customizing aesthetics
```r
ggplot(mtcars, aes(x = mpg, y = disp, color = cyl)) +
  geom_point()
```


Anything you put in the ggplot() function can be seen by any geom layers that you add (i.e., these are universal plot settings). This includes the x- and y-axis you set up in aes().

You can also specify aesthetics for a given geom independently of the aesthetics defined globally in the ggplot() function.

The + sign used to add layers must be placed at the end of each line containing a layer. If, instead, the + sign is added in the line before the other layer, ggplot2 will not add the new layer and will return an error message.






## Customization

`ggplot2` allows for extensive customization. For example, you can add a title to your plot and label your axes as follows:

```r
ggplot(data = mtcars, aes(x = mpg, y = disp)) +
  geom_point() +
  labs(title = "Miles Per Gallon vs. Displacement",
       x = "Miles Per Gallon",
       y = "Displacement")
```

## Conclusion

`ggplot2` is a versatile tool for creating high-quality visualizations in R. With its simple syntax and extensive customization options, it's a valuable addition to any data scientist's toolbox.

---

Remember to replace the code snippets with your own examples and explanations. Happy coding! 🚀.

Source: Conversation with Bing, 22/01/2024
(1) en.wikipedia.org. https://en.wikipedia.org/wiki/Ggplot2.
