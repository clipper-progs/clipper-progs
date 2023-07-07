# clipper-progs
This repository contains the suite of clipper-based programs used within CCP4.

## Project Structure

```mermaid
graph TD
A[clipper-progs] --> B(buccaneer)
A[clipper-progs] --> C(ctruncate)
A[clipper-progs] --> D(nautilus)
A[clipper-progs] --> E(parrot)
A[clipper-progs] --> F(pirate)
A[clipper-progs] --> G(sheetbend)
A[clipper-progs] --> H(sloop)
```

## Building
### First Time
    git clone https://github.com/clipper-progs/clipper-progs.git
    cd clipper-progs
    git submodule update --init --recursive

### To get the latest changes
    git pull
    git submodule update --recursive
