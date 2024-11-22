# Advanced Data Visualization with ggplot2

**ggplot2** is a powerful and flexible R package for creating stunning and informative visualizations. It implements the grammar of graphics, making it easy to understand and implement various data visualizations. This expanded tutorial is designed for a 2-hour class, providing a comprehensive guide to **ggplot2** suitable for self-learning. We will cover both basic and advanced features, providing detailed explanations for each `geom` and `stat` used.

This tutorial assumes a working knowledge of **R** and **RStudio**.

### References
- [ggplot2 Reference](https://ggplot2.tidyverse.org/reference/)
- [ggplot2 Specs](https://ggplot2.tidyverse.org/articles/ggplot2-specs.html)
- [R Graphics Cookbook](https://r-graphics.org/)

## Installation
To begin using **ggplot2**, install it from **CRAN** using:

```r
install.packages("ggplot2")
```

Alternatively, you can use **ggplot2** from the **tidyverse** package:

```r
install.packages("tidyverse")
```

## Basic Usage
First, load **ggplot2** into your **R** environment:

```r
library(ggplot2)
```

The main function is `ggplot()`, which serves as a base layer for adding additional components to your plot.

```r
ggplot(data = <DATA>, mapping = aes(<MAPPINGS>)) + <GEOM_FUNCTION>()
```

- **Data**: The raw data to be plotted.
- **Aesthetics (`aes()`)**: Define attributes like position, color, size, shape, etc.
- **Geometries (`geom_`)**: Specify the shapes that represent your data, like points or lines.

## Overview of the `mtcars` Dataset
In this tutorial, we use the **mtcars** dataset. To preview it:

```r
head(mtcars)
```

The **mtcars** dataset contains information about fuel consumption and other aspects of automobile design and performance for 32 automobiles (1973–74 models). The dataset includes the following variables:

- **mpg**: Miles per gallon
- **cyl**: Number of cylinders
- **disp**: Displacement (in cubic inches)
- **hp**: Gross horsepower
- **drat**: Rear axle ratio
- **wt**: Weight (in 1000 lbs)
- **qsec**: 1/4 mile time
- **vs**: Engine shape (0 = V-shaped, 1 = straight)
- **am**: Transmission (0 = automatic, 1 = manual)
- **gear**: Number of forward gears
- **carb**: Number of carburetors

## Basic Plot with `geom_point`
The simplest plot you can create is a scatter plot, which shows the relationship between two variables:

```r
ggplot(data = mtcars, mapping = aes(x = wt, y = mpg)) +
  geom_point()
```

- **`geom_point()`**: Adds a layer of points to the plot. It is useful for showing the relationship between two numeric variables.
  - **`x` and `y`**: These are the variables that determine the positions of the points on the x-axis and y-axis respectively.

In this example, **`wt`** (weight of the car) is plotted on the x-axis and **`mpg`** (miles per gallon) is plotted on the y-axis, helping us see how car weight impacts fuel efficiency.

### Customizing `geom_point`
You can customize the appearance of points by using additional aesthetics within `aes()` such as color, size, and shape:

```r
ggplot(data = mtcars, mapping = aes(x = wt, y = mpg, color = factor(cyl), size = hp)) +
  geom_point()
```

- **`color`**: Differentiates points by a categorical variable (e.g., **cyl** for the number of cylinders). Here, we use `factor(cyl)` to ensure it is treated as a categorical variable.
- **`size`**: Modifies the point size based on a continuous variable (e.g., **hp** for horsepower).

In this example, cars with different numbers of cylinders are represented by different colors, and their horsepower determines the size of the points, providing more information at a glance.

### Adding Layers
Layers can be added to enhance the plot, such as adding trend lines or statistical summaries:

```r
ggplot(data = mtcars, mapping = aes(x = wt, y = mpg)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE)
```

- **`geom_smooth()`**: Adds a smoothing line to represent trends in the data.
  - **`method = "lm"`**: Specifies that a linear model should be used to fit the line, helping to visualize the trend.
  - **`se`**: Controls whether to display the confidence interval around the smoothing line. Setting `se = FALSE` removes the shaded confidence region.

This combination of layers helps to identify overall trends in the data, with **`geom_point()`** showing individual data points and **`geom_smooth()`** providing a fitted trend line.

## Advanced Geoms

### `geom_line()`
Used for line charts, typically to show trends over time.

```r
ggplot(data = economics, aes(x = date, y = unemploy)) +
  geom_line(color = "blue")
```

- **`geom_line()`**: Plots lines, useful for time-series data.

### `geom_bar()`
Used for creating bar plots. There are two main types: **count-based** and **identity-based**.

#### Count-Based Bar Plot

```r
ggplot(mtcars, aes(x = factor(cyl))) +
  geom_bar()
```

- **`geom_bar()`**: By default, counts the occurrences of each category.

#### Identity-Based Bar Plot
If you have pre-counted data, use `stat = "identity"`:

```r
# Assume df is a dataset with columns category and count
ggplot(df, aes(x = category, y = count)) +
  geom_bar(stat = "identity", fill = "skyblue")
```

- **`stat = "identity"`**: Directly uses values from the dataset.
- **`fill`**: Specifies the fill color of bars.

## Faceting
Faceting is a powerful way to split your data into multiple panels based on categorical variables, allowing for easy comparison between different groups.

```r
ggplot(mtcars, aes(x = mpg, y = disp)) +
  geom_point() +
  facet_wrap(~cyl)
```

- **`facet_wrap()`**: Creates sub-plots for each unique value of a variable (e.g., **cyl** for number of cylinders). This allows you to see how the relationship between **mpg** and **disp** changes across different cylinder configurations.

### `facet_grid()`
To create a matrix of plots based on two variables, you can use `facet_grid()`:

```r
ggplot(mtcars, aes(x = wt, y = mpg)) +
  geom_point() +
  facet_grid(gear ~ cyl)
```

- **`facet_grid()`**: Arranges plots in a grid format based on two variables, allowing for detailed comparative analysis. Here, rows represent different values of **gear** and columns represent different values of **cyl**.

## Statistical Transformations (`stat_` Functions)
Statistical transformations help represent your data in different ways, such as histograms or summary statistics.

### `geom_histogram()`
Histograms are useful for displaying the distribution of a single numeric variable, giving insight into the frequency of values within different ranges.

```r
ggplot(mtcars, aes(x = mpg)) +
  geom_histogram(binwidth = 2, fill = "orange", color = "black")
```

- **`geom_histogram()`**: Divides data into bins and shows the frequency.
- **`binwidth`**: Adjusts the width of bins, controlling granularity.

### `geom_boxplot()`
Useful for showing the distribution, median, and potential outliers of numeric data.

```r
ggplot(mtcars, aes(x = factor(cyl), y = mpg)) +
  geom_boxplot(notch = TRUE, fill = "lightblue")
```

- **`geom_boxplot()`**: Displays a box plot, summarizing data distribution.
- **`notch = TRUE`**: Adds notches, indicating the confidence interval around the median.

## Adding Annotations and Labels
Annotations are key to making plots informative. You can add titles, labels, and text annotations to your plot.

### Adding Titles and Axis Labels

```r
ggplot(mtcars, aes(x = mpg, y = disp)) +
  geom_point() +
  labs(title = "Miles Per Gallon vs. Displacement",
       x = "Miles Per Gallon",
       y = "Displacement")
```

### Adding Text Labels to Points
To make your plots more informative, you can add labels to each data point:

```r
ggplot(mtcars, aes(x = mpg, y = disp)) +
  geom_point() +
  geom_text(aes(label = rownames(mtcars)), vjust = 1, hjust = 1)
```

- **`geom_text()`**: Adds text labels to each point, which helps in identifying specific data points directly on the plot.
- **`vjust`** and **`hjust`**: Control the vertical and horizontal positioning of the text labels relative to the points.

## Customizing Themes
Themes control the overall appearance of your plot, such as background color, grid lines, and fonts.

### Applying Predefined Themes

```r
ggplot(mtcars, aes(x = mpg, y = disp)) +
  geom_point() +
  theme_minimal()
```

- **`theme_minimal()`**: Applies a clean, minimalistic style.

### Creating Custom Themes

```r
ggplot(mtcars, aes(x = mpg, y = disp)) +
  geom_point() +
  theme(panel.background = element_rect(fill = "white", color = "black"),
        panel.grid.major = element_line(color = "gray"))
```

- **`element_rect()`**: Defines a rectangular element, used here for the plot background.
- **`element_line()`**: Customizes grid lines.

## Interactive Visualizations with `plotly`
You can convert static **ggplot2** plots to interactive ones using the **plotly** package.

```r
library(plotly)
g <- ggplot(mtcars, aes(x = wt, y = mpg)) +
  geom_point()
ggplotly(g)
```

- **`ggplotly()`**: Converts a `ggplot` object to an interactive plot.

## Practice Exercise
Create a multi-faceted plot of the **mpg** variable against **wt**, with each facet representing a different **gear** value from the **mtcars** dataset. Customize the color by **cyl** and add a linear trend line. Use a suitable theme to make the plot publication-ready.

### Solution

```r
ggplot(mtcars, aes(x = wt, y = mpg, color = factor(cyl))) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~ gear) +
  theme_bw() +
  labs(title = "Weight vs. MPG by Gear Number",
       x = "Weight (1000 lbs)",
       y = "Miles Per Gallon",
       color = "Number of Cylinders")
```

## Summary
In this expanded tutorial, we covered:
- Basic and advanced plotting with **ggplot2**
- How to add layers and statistical transformations
- Customizing aesthetics, scales, and themes
- Creating faceted plots for better data comparison
- Adding informative annotations and labels

The power of **ggplot2** lies in its flexibility and its layered approach to building plots. With the knowledge you've gained, you should be well-equipped to create publication-quality visualizations that communicate your data clearly and effectively.

Continue practicing by applying these techniques to your own datasets and exploring the extensive customization options available in **ggplot2**.

