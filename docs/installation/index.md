<br>

# Installing G4X-helpers

These pages explain how you can obtain G4X-helpers and set up your system to use it. 
Installation can be done either via source installation or utilizing a Docker image.  
Which route to choose depends on your use-case.

<br>

### 1. source installation
---

**Please head to the [source installation](./source.md) section if you want to:**

+ use the `G4X-helpers` CLI directly on your machine
+ use the package beyond the CLI features exposed in the Docker image
+ develop or debug the codebase
+ integrate pieces of the library into your own Python environment

### 2. Docker setup
---

**Please head to the [Docker setup](./docker.md) section if you want to:**

+ avoid creating or managing local Python environments
+ guarantee reproducibility (everyone runs the exact same image)
+ run G4X-helpers easily on different systems
+ stay isolated and want to ensure that nothing leaks to or depends on your system
 
Things to know:  

+ you still need to install Docker (or Apptainer, Podman ... ) if you haven't already.
+ first pull can be large; subsequent runs are fast thanks to caching layers.

--8<-- "_partials/end_cap.md"
