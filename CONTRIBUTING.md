# Contributing to ParaDiGM

First of all, thank you for considering contributing to ParaDiGM! 

## Development Model
Please note that the core development of ParaDiGM takes place on an internal GitLab server at **ONERA**. This GitHub repository serves as the official mirror and distribution point for the community.

## How to Contribute

### 1. Bug Fixes and Minor Improvements
For small contributions, such as fixing a typo, updating documentation, or correcting a minor bug:
* You can submit a **Pull Request** directly on this GitHub repository.
* Ensure your code follows the existing coding standards (see below).
* Your **Pull Request** will be reviewed by the core maintainers and, if accepted, integrated into the internal repository before being pushed back to GitHub.

### 2. Major Contributions
For significant features, algorithmic evolutions or architectural changes:
* **Requirement:** Please **open an Issue** on GitHub first to discuss your proposal with the development team. 
* This ensures that your work aligns with the project roadmap and avoids duplicate efforts.

## Coding Standards

To maintain a clean and readable codebase, ParaDiGM follows professional standards. While we do not enforce a strict internal manual, our C/C++ code closely adheres to the [Google C++ Style Guide](https://google.github.io/styleguide/cppguide.html).

### Specific Rules:
* **Naming (C++):** Use `camelCase` for functions and variables.
* **Naming (C):** Use `snake_case` for functions and variables to maintain consistency with standard C libraries.
* **Indentation:** We use **2 spaces** per indentation level (no tabs).
* **Braces:** The opening brace `{` must be placed **at the end of the line** (not on a new line), as follows:
  ```cpp
  if (condition) {
    // code
  }
  
## Attribution
By contributing, you agree that your contributions will be licensed under the project's [LICENSE](LICENSE) (LGPL). All contributors will be acknowledged in the documentation.