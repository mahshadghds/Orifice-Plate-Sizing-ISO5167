# Orifice Plate Sizing Tool (BS EN ISO 5167)

A computational engineering tool developed in **C++** implementing the **BS EN ISO 5167** standard for differential pressure flowmeter sizing and fluid flow calculation.

## 📄 Documentation & Technical Report
* **[Read the Full Technical Validation Report (PDF)](Validation_Report.pdf)** — Comprehensive mathematical derivation, iterative sizing methodology, and validation against industrial datasheets (achieving less than **0.1%** deviation).
* **[View High-Resolution Algorithmic Flowchart (PDF)](Algorithm_Flowchart.pdf)**

---

## 🔄 Algorithmic Logic & Flowchart
The iterative calculation procedure for computing the discharge coefficient ($C$), expansion factor ($\epsilon$), Reynolds number ($Re_D$), and orifice bore diameter ($d$) is illustrated below:

![Algorithm Flowchart](flowchart.png)

---

## 💻 Source Code & Usage
The full implementation is located in [`orifice_sizing.cpp`](orifice_sizing.cpp).

To compile and run the program using `g++`:
```bash
g++ -O3 orifice_sizing.cpp -o orifice_sizing
./orifice_sizing
