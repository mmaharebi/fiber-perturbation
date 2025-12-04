# Code Boxes Guide

The template provides beautiful, syntax-highlighted code boxes that match the design system.

## Available Code Environments

### 1. **Generic Code Box with Title**

```latex
\begin{codebox}[language=Python]{Custom Title}
def hello_world():
    print("Hello, World!")
\end{codebox}
```

**Features:**
- Custom title
- Line numbers
- Syntax highlighting
- Breakable across pages
- Matches design system colors

### 2. **Inline Code Block (No Title)**

```latex
\begin{codeblock}[language=Python]
x = 42
print(x)
\end{codeblock}
```

**Features:**
- No title, minimal styling
- No line numbers
- Good for short snippets
- Lighter frame color

### 3. **Python Code Box**

```latex
\begin{pythoncode}
import numpy as np

def calculate_mean(data):
    """Calculate the mean of a dataset."""
    return np.mean(data)

# Example usage
data = [1, 2, 3, 4, 5]
result = calculate_mean(data)
print(f"Mean: {result}")
\end{pythoncode}
```

### 4. **Java Code Box**

```latex
\begin{javacode}
public class HelloWorld {
    public static void main(String[] args) {
        System.out.println("Hello, World!");
    }
}
\end{javacode}
```

### 5. **C++ Code Box**

```latex
\begin{cppcode}
#include <iostream>

int main() {
    std::cout << "Hello, World!" << std::endl;
    return 0;
}
\end{cppcode}
```

## Supported Languages

The `listings` package supports many languages:

- **Programming**: Python, Java, C, C++, C#, JavaScript, TypeScript, PHP, Ruby, Go, Rust, Swift, Kotlin
- **Scripting**: Bash, PowerShell, Perl, Lua
- **Markup**: HTML, XML, LaTeX
- **Database**: SQL, MySQL, PostgreSQL
- **Math/Scientific**: MATLAB, R, Mathematica
- **Functional**: Haskell, Lisp, Scheme, Erlang, OCaml
- **Assembly**: x86, ARM

## Customization Options

### Change Language

```latex
\begin{codebox}[language=Java]{Matrix Multiplication}
// Java code here
\end{codebox}
```

### Remove Line Numbers

```latex
\begin{codebox}[language=Python, numbers=none]{Without Line Numbers}
# Python code here
\end{codebox}
```

### Change Starting Line Number

```latex
\begin{codebox}[language=Python, firstnumber=10]{Starts at Line 10}
# This will be line 10
# This will be line 11
\end{codebox}
```

### Emphasize Specific Lines

```latex
\begin{codebox}[language=Python, emphstyle=\color{warning-red}]{Highlighted Lines}
def example():
    x = 1  # Normal line
    y = 2  # Normal line
\end{codebox}
```

### Custom Tab Size

```latex
\begin{codeblock}[language=Python, tabsize=2]
def foo():
  return 42
\end{codeblock}
```

## Color Scheme

The code boxes use your design system colors:

- **Background**: `bg-light` (off-white)
- **Frame**: `example-slate` (slate gray) for titled boxes
- **Frame**: `divider` (light gray) for inline blocks
- **Keywords**: `primary-main` (royal blue), bold
- **Comments**: `text-secondary` (gray), italic
- **Strings**: `success-green` (green)
- **Line numbers**: `text-secondary` (gray), tiny

## Advanced: Using Minted (Optional)

For even better syntax highlighting, you can use `minted` instead:

### 1. Install Pygments

```bash
pip install Pygments
```

### 2. Uncomment in `config/packages.tex`

```latex
\usepackage{minted}
```

### 3. Compile with `-shell-escape`

```bash
pdflatex -shell-escape main.tex
```

### 4. Use minted environments

```latex
\begin{minted}{python}
def hello():
    print("Hello with Pygments!")
\end{minted}
```

## Tips

1. **For short snippets**: Use `\verb|code|` or `\texttt{code}` inline
2. **For terminal output**: Use `codeblock` without syntax highlighting
3. **For algorithms**: Consider using the `algorithm2e` package instead
4. **For pseudocode**: Use `algorithmicx` package

## Examples in the Template

See `sections/04_advanced_features.tex` for live examples of code boxes.

## Troubleshooting

**Problem**: Code doesn't fit in box  
**Solution**: Add `breaklines=true` (already default)

**Problem**: Strange characters in code  
**Solution**: Use `\verb` or escape special LaTeX characters

**Problem**: Need Unicode support  
**Solution**: Switch to XeLaTeX/LuaLaTeX or use `minted`

---

**Note**: Code boxes are breakable across pages by default. To prevent breaking, add `breakable=false` to the options.
