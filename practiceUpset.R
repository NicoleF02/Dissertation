library(ggplot2)
library(ComplexUpset)

movies = as.data.frame(ggplot2movies::movies)

genres = colnames(movies)[18:24]
movies[genres] = movies[genres] == 1

movies[movies$mpaa == '', 'mpaa'] = NA
movies = na.omit(movies)


upset(movies, genres, name='genre', width_ratio=0.1)
