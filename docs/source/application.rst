PyThea Application
==================

The Graphical User Interface
----------------------------

Performing a 3D reconstruction is mainly an interactive process. To reconstruct an event, the user aims to achieve the best fit of a geometrical model to multi-viewpoint coronal observations by adjusting a set of geometrical parameters. However, without a graphical user interface (GUI), this process becomes nearly impossible. PyThea provides a modern solution for conducting comprehensive analyses of events. Its user-friendly GUI, built on Streamlit, facilitates efficient interaction with the reconstruction process.

The PyThea web application comprises two main pages as shown in :numref:`main-app-page`. The main page (left panel) appears when running the application for the first time. On this page, you can initiate the fitting process by selecting the date of the event and the fitting model to be used. Then the application initiates and a second page is loaded with the fitting sliders and the imaging views. Each main page comprises two primary vertical panels as shown in :numref:`main-app-page`, namely Panel 1 and 2. The left panel serves as a hub for your input widgets, allowing you to interactively change the parameters and control the reconstruction process. On the other hand, the right panel is dedicated to displaying data elements, providing a visual feedback on the reconstruction results. This separation of input and output elements enhances the usability and clarity of the application interface, enabling you to navigate through the reconstruction process with ease. In essence, PyThea revolutionizes the 3D reconstruction workflow by combining advanced algorithms with intuitive interface design, making complex analyses accessible to a broader audience.

.. figure:: ./images/main_page.png
   :name: main-app-page
   :width: 800px
   :align: center

   Two views of PyThea‘s web application. The left panel shows the starting page of the application and the right panel shows the main fitting page.

Initializing the Fitting Process
--------------------------------

When the application runs for the first time, the starting page appears and you have to select the date and event to process and the geometrical model to use. A view of the starting page is shown in :numref:`main-app-page`. For the selection of the date and event selection, three options are available to initiate the fitting process. With the 'Event' mode you can select a date and a solar flare event associated with the event that you want to perform the geometrical fitting. If the CME or shock do not have an associated flare then there is a second 'Manual' option. With this option, you can select the date and the event (CME, shock, flare). A third option is to provide a fitting file and the application will initiate with the fitted event and the application will load the previous fitting parameters. Next, you will select the geometrical model for the fitting process. After this selection the application will search and download available data and the main geometrical fitting page will appear.

Main Geometrical Fitting Page
-----------------------------

When the main geometrical fitting page appears you can start the fitting process. In this page you can adjust the parameters of the geometrical model, change the images and the selected imager, and select among various visualization options. An example of the main geometrical fitting page is shown in :numref:`main-app-page` (right panel).

.. figure:: ./images/app_slider_details.png
   :name: app_slider_details
   :width: 800px
   :align: center

   These panels show an example of the input widgets that can be used to provide input parameters to the application.

Widgets Panel
~~~~~~~~~~~~~

You can easily customize your visualization using various input widgets found on the left panel of this page. From sliders to radio buttons and drop-down menus, these widgets offer a range of options for visualizing your data, adjusting fitted parameters, selecting specific images, and performing image processing tasks. An example of the different options is shown in :numref:`app_slider_details`.

In the first container (panel a), you'll find options related to the coordinate system and axis representation for fitting parameters, and visualization options for the model, as well as, visualization of the fitting table and kinematics. Switching between coordinate systems allows you to specify model coordinates in Carrington or Stonyhurst systems, offering flexibility in your analysis. Changing the coordinate system of the geometrical model you can provide the longitude and latitude of the model in Carrington or Stonyhurst coordinates. Similarly, the axis representation provides alternative ways to input geometrical parameters for the model. You can also select whether to display the full mesh or just the skeleton of the geometrical model within the images. Additionally, the 'View Fitting Table' and 'View Kinematics Plots' options enable you to visualize the fitting table containing the geometrical parameters and kinematic plots resulting from the fittings, respectively. These features provide comprehensive insights into your fittings and the kinematics of the model during the analysis.

The next container has the sliders with the geometrical parameters. Adjusting the parameter values of the geometrical model will update its view at the selected imager to the new state. Adjust the parameters until there is a good fit of the model to the CME or shock wave that appears on the images. When the fitting is ready you can store it to the fitting table by pressing the 'Store Fit' button. Storing a new fit to the table will update the table and the kinematic plots. If the fitting for the selected time and image already exists in the table then the values will be updated.

The 'Imaging menu' (panel c) provides options to select more imagers and download and load the imaging data. With the 'Time Range' slider you can extend the time interval of the image loading. The preselected time interval is one hour before and after the flare maximum of the time selected with the manual mode. You can also select among three different image processing options in panel c, namely running and base difference images, and plain images. There is also an option to clip the geometrical model on the image limits. Additionally, the option of 'Supplementary Imaging' visualizes two near-simultaneous images from other imagers. This can be used to perform triangulation and tight constrain the geometrical model using multi-viewpoint images. The image's climits slider provide also the option to change the colormap limits of the images shown.

Imaging and Results Panel
~~~~~~~~~~~~~~~~~~~~~~~~~

The right panel of the main geometrical fitting page is the imaging panel. Here appear the images for the selected imager with the geometrical model overlaid onto the images. You can select the imager to visualize the available data from a drop drop-down menu and also use a slider to select the different times of the images. Additionally, if the 'Supplementary Imaging' has been selected two nearly co-temporal images from other instruments will appear. You can use the slider on the bottom of the two images to select among the different loaded imagers. In this panel, the fitting table and the kinematic plots can appear on the bottom of the images, if these options have been selected. We give more details in the following Section.

Fitting Table and Kinematics Plot
---------------------------------
