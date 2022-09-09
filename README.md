# Plant Generator

This extension lets you generate plants procedurally. There are a load of settings to change, so have fun :)!

How to use:
1. (optional) Edit settings to change the shape of the plant
2. (optional) Change the "Plant Group" field. This will put generated plants under this folder
3. Click the "Draw Plant" button
4. Click around in the scene to place plants on meshes (note - does not work on ground plane objects because mesh raycast doesn't hit those)
5. Esc stops plant draw mode

If you mess up a settings, settings can be reset to defaults with the "Reset" buttons

Plant settings can be changed in real time.
Simply select your plants (or a parent of plants) and edit the settings. Plants will automatically update.


Note - the performance is poor as of version 1.0.0. This is because it's written in python and I had no time to optimize it.