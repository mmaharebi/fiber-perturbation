# Design System Documentation

This template uses a comprehensive, professional design system that ensures visual consistency across all elements.

## 🎨 Color Palette

### Academic Professional Theme

The template uses a carefully curated color palette based on Material Design principles, optimized for academic and technical documents.

#### Primary Colors (Main Content)
- **Primary Dark** `RGB(25, 42, 86)` - Deep Navy for headings
- **Primary Main** `RGB(41, 98, 255)` - Royal Blue for links and accents
- **Primary Light** `RGB(100, 149, 237)` - Cornflower Blue for lighter accents

#### Secondary Colors (Supporting Elements)
- **Secondary Dark** `RGB(0, 77, 64)` - Teal Dark
- **Secondary Main** `RGB(0, 150, 136)` - Teal
- **Secondary Light** `RGB(77, 182, 172)` - Teal Light

#### Semantic Colors (Theorem Boxes)
- **Theorem Blue** `RGB(25, 118, 210)` - For theorems, lemmas, propositions, corollaries
- **Definition Teal** `RGB(0, 137, 123)` - For definitions and foundational concepts
- **Example Slate** `RGB(84, 110, 122)` - For examples and illustrations
- **Remark Purple** `RGB(123, 31, 162)` - For remarks, notes, and observations
- **Question Amber** `RGB(255, 143, 0)` - For questions and problems
- **Success Green** `RGB(56, 142, 60)` - For solutions and answers
- **Warning Red** `RGB(211, 47, 47)` - For warnings and critical alerts
- **Important Orange** `RGB(245, 124, 0)` - For important information

#### Neutral Colors
- **Text Primary** `RGB(33, 33, 33)` - Main body text
- **Text Secondary** `RGB(117, 117, 117)` - Secondary text, captions
- **Background Paper** `RGB(255, 255, 255)` - White background
- **Background Light** `RGB(250, 250, 250)` - Off-white for code blocks
- **Divider** `RGB(224, 224, 224)` - Light gray for rules and borders

### Using Colors in Your Document

Colors are automatically applied to:
- Section headings
- Theorem boxes
- Links and references
- Lists
- Code blocks
- Tables

To use colors manually:
```latex
{\color{primary-main}This text is royal blue}
{\color{theorem-blue}Theorem-colored text}
```

## 📝 Typography

### Font Configuration

**Default**: Latin Modern (enhanced Computer Modern)

**Optional Professional Fonts** (uncomment in `config/theme.tex`):
- `tgtermes` - Times-like serif
- `tgpagella` - Palatino-like serif  
- `libertine` - Linux Libertine serif
- `FiraSans` - Modern sans-serif

### Heading Styles

- **Part**: Huge, bold, deep navy with spacing
- **Section**: Large, bold, deep navy with underline
- **Subsection**: Large, bold, royal blue
- **Subsubsection**: Normal size, bold, teal

### Text Spacing

- **Line spacing**: 1.15 (slightly increased for readability)
- **Paragraph indent**: 0pt (no indentation)
- **Paragraph skip**: 0.8em (space between paragraphs)

## 📦 Theorem Boxes

### Numbered Environments

All theorem environments use consistent styling with semantic colors:

```latex
\begin{theorem}{Title}{label}
  Content here
\end{theorem}

\begin{definition}{Title}{label}
  Content here
\end{definition}

\begin{lemma}{Title}{label}
  Content here
\end{lemma}

\begin{example}{Title}{label}
  Content here
\end{example}

\begin{remark}{Title}{label}
  Content here
\end{remark}
```

### Utility Boxes (Non-numbered)

For special purpose content without numbering:

```latex
\begin{warningbox}
  Critical information
\end{warningbox}

\begin{importantbox}
  Key points to remember
\end{importantbox}

\begin{questionbox}
  What is the derivative of x²?
\end{questionbox}

\begin{answerbox}
  The answer is 2x
\end{answerbox}

\begin{solutionbox}
  Detailed solution steps...
\end{solutionbox}

\begin{tipbox}
  Helpful hint
\end{tipbox}
```

### Box Styling Parameters

All boxes share consistent styling:
- **Box rule**: 0.8pt
- **Corner radius**: 3pt
- **Padding**: 8pt on all sides
- **Colored frame** with 5% tinted background
- **White title text** on colored background

## 📋 Lists

### Itemized Lists

- First level: Royal blue bullets (•)
- Second level: Teal circles (○)

### Enumerated Lists  

- First level: Royal blue numbers (1., 2., 3.)
- Second level: Teal letters (a), b), c))

### Custom Spacing

- Item separation: 0.3em (first level)
- Item separation: 0.2em (second level)
- Parse separation: 0em

## 🔗 Links and References

- **Internal links**: Royal blue (`primary-main`)
- **URLs**: Royal blue (`primary-main`)
- **Citations**: Teal (`secondary-main`)
- **File links**: Royal blue (`primary-main`)

All links have no borders and break across lines when needed.

## 📊 Tables

- **Rule color**: Light gray divider
- **Row spacing**: 1.3× (increased for readability)
- **Caption labels**: Bold, deep navy
- **Caption text**: Primary text color

## 💻 Code Listings

- **Background**: Off-white (`bg-light`)
- **Frame**: Light gray divider
- **Basic text**: Primary text color, small monospace
- **Comments**: Secondary text, italic
- **Keywords**: Royal blue, bold
- **Strings**: Success green
- **Line numbers**: Tiny, secondary text

## 🖼️ Figures and Captions

- **Caption labels**: Bold, deep navy
- **Caption text**: Primary text color, small
- **Subcaption labels**: Bold, teal
- **Subcaption text**: Secondary text, footnote size

## 📄 Table of Contents

- **Section entries**: Bold, deep navy
- **Subsection entries**: Primary text color
- **Page numbers**: Matching text colors
- **Leaders**: Dotted lines

## 🎯 Custom Commands

The theme provides helpful custom commands:

```latex
% Highlighted important text
\highlight{This is important!}

% Subtle secondary text
\subtle{Additional information}

% Custom horizontal rule
\customhrule
```

## 🎨 Customization

### Changing the Color Scheme

Edit `config/theme.tex` to modify colors:

1. **For a warmer palette**: Adjust primary colors to browns/oranges
2. **For a cooler palette**: Use blues/purples/teals
3. **For grayscale**: Reduce saturation of all semantic colors

### Changing Fonts

In `config/theme.tex`, uncomment one of the font packages:

```latex
\usepackage{tgpagella}  % Palatino-like
```

### Enabling Headers/Footers

Uncomment the `fancyhdr` section in `config/theme.tex`:

```latex
\usepackage{fancyhdr}
\pagestyle{fancy}
% ... (styling code)
```

### Adjusting Box Styles

Modify the `\boxstyle` command in `config/theorems.tex`:

```latex
\newcommand{\boxstyle}{
  fonttitle=\bfseries,
  boxrule=1.0pt,        % Thicker border
  arc=5pt,              % More rounded corners
  left=10pt,            % More padding
  % ...
}
```

## 🔧 Technical Notes

### Package Loading Order

The correct order is critical:

1. **packages.tex** - Core LaTeX packages
2. **theme.tex** - Colors and styling (must come after xcolor, hyperref)
3. **macros.tex** - Custom commands
4. **theorems.tex** - Theorem boxes (uses theme colors)
5. **metadata.tex** - Document info

### Color Dependencies

- `xcolor` must be loaded before `theme.tex`
- `hyperref` must be loaded before `theme.tex` (for link styling)
- `tcolorbox` must be loaded before `theorems.tex`

### Compatibility

This design system is compatible with:
- LaTeX, pdfLaTeX, XeLaTeX, LuaLaTeX
- Standard document classes (article, report, book)
- Common packages (amsmath, tikz, biblatex, etc.)

## 📚 Best Practices

1. **Consistency**: Use the predefined theorem boxes rather than custom tcolorboxes
2. **Accessibility**: The color palette has good contrast ratios for readability
3. **Printing**: Colors are chosen to look good in both digital and printed form
4. **Semantic meaning**: Use colors according to their semantic purpose

## 🚀 Quick Start

To use the design system in a new document:

```latex
\documentclass{article}

\input{config/packages}
\input{config/theme}
\input{config/macros}
\input{config/theorems}
\input{config/metadata}

\begin{document}
% Your content here
\end{document}
```

## 📖 Examples

See the demonstration sections in the template for examples of:
- `sections/03_theorem_boxes.tex` - All theorem box styles
- `sections/04_advanced_features.tex` - Advanced formatting
- `sections/05_complete_example.tex` - Complete document examples

---

**Version**: 1.0  
**Last Updated**: November 2025  
**Theme**: Academic Professional
