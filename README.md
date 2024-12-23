
# MATH521 - Numerical Analysis of PDE

This git repository contains some course materials for `MATH521 - Numerical Analysis of PDE`, taught at UBC, by C Ortner. I will update/change the content throughout term and announce updates in class. Feel free to try it out but be aware that it might change throughout term. 

### Get started 

- All materials are provided as Julia scripts. To install Julia, follow the instructions at [https://julialang.org](https://julialang.org). I suggest to use at least Julia 1.9, but better 1.10. 
- Clone this `git` repository to your local computer, say to the folder `.../math521`. If you don't know how to use git, Atlassian makes [some great tutorials](https://www.atlassian.com/git/tutorials). If you don't want to learn git, then you can also just download the repository, but it will make it awkward to update it.
- You can setup your environment by starting Julia in the folder `.../math521`, then start the Julia REPL, type `]` to switch to the package manager and then 
``` 
activate .
up
```
Hopefully this should get you fully set up. 
- Now press backspace to get back to the main REPL mode and type 
```julia 
using IJulia 
IJulia.notebook() 
```
which should start a jupyter server and from there you can open the Jupyter notebooks. 

### Notes 

- If haven't used Julia in the past, read the [Get Started](https://julialang.org/learning/getting-started/) pages. There are also plenty of excellent tutorials as at [`.../learning`](https://julialang.org/learning/).
- It is strictly speaking not required to learn Julia for this course. You can equally use Python (possibly also Matlab, but I'm less sure about that ...). If you port some of the course material to another language, let me know and maybe we can add it to the course pages. 
- If the instructions in this readme are insufficient, confusing or incorrect, please file an issue or come to my office hour or make a PR to fix them. 



