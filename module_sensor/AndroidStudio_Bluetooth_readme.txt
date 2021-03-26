Prerequisites
To implement the project, the following software must be installed:

Android Studio [https://developer.android.com/distribute]
Arduino [https://www.arduino.cc/en/software]
Unity [https://store.unity.com/download]

For this project we will be using the following builds and versions:

Android Studio
* Gradle Build: 3.5.1
* Gradle Version: 6.1.1

Arduino
* ESP32 Add-on in Arduino IDE

Unity
* Version 2018.4.24f1

The code for the project can be cloned from the Github repository (https://github.com/mylespedronan/prosthetic)
Building the Project


1. Open Android Studio and open the Prosthetic project


2. Check that both the Gradle Build and Gradle Version match the prerequisites above by checking Prosthetic > build.gradle (Project: Prosthetic) and Prosthetic > gradle > wrapper > gradle-wrapper.properties


Build.gradle file


gradle-wrapper.properties


3. Open the Arduino IDE and install the ESP32 add on by following the tutorial on (https://randomnerdtutorials.com/installing-the-esp32-board-in-arduino-ide-windows-instructions/)


4. Connect Arduino to the correct COM port in Tools and compile the code onto the ESP32 board. For this example we will be using the handApp.ino file.


Connecting COM port


Compile button


a. If the terminal shows “Connecting..._____”, hold down the Boot button on the ESP32 until the terminal starts to compile the code. If the code is successfully built, the terminal will show “Leaving… Hard resetting via RTS pin…”. 


Compiling before pressing the boot button on the ESP32 board


Compiling after the boot button


b. You can verify the code has been properly compiled by opening the serial monitor on the Arduino IDE and checking if the program is working as intended. The baud rate used is 115200. With the serial monitor opened, press the reset button on the ESP32 and verify with the monitor that the device can be paired with bluetooth. If successful, you can close the serial monitor.


Arduino Serial Monitor


5. Open the Unity project and change the window to view the console


Hand Model in Unity


6. Press the play button in the top middle of the Unity window and verify there are no errors with the code. If there are no errors, the console window will be empty 


Play button pressed


7. Return to Android Studio and connect your phone via USB. Run the app by pressing the play button in the toolbar.


Run the application on your phone


8. Once the application compiles you will be brought to the homepage of the app


Home screen


9. Using the side bar on the upper left hand of the screen we can bring up the other available layouts. Here we will click on “Connect Device”


Side Bar


10. On the Connect Device page, we will need to connect to the ESP32 in Bluetooth settings before we can connect using the application. To do this, click on the button on the top right of the screen to bring up “Bluetooth Settings”


Bluetooth Settings


11. In the Bluetooth Settings, select the ESP32 board and once connected, go back to the app. If successful, the ESP32 board will now be selectable in the Connect Devices page


Connecting to ESP32 in Bluetooth Settings


12. Click on the ESP32 board to go to the terminal. Here the terminal will attempt to connect to the ESP32 and upon successful completion, will verify that it is connected


Connect to ESP32 in app and Terminal Window


13. Using the voice button in the terminal, you will be able to speak any of the commands below. Allow the application to record audio.


Allowing app to record audio


List of Commands for Hand Model


14. For this example we will use “Activate large diameter”. Press the play button and say, “Activate large diameter”. The application will repeat the command said and show the character array being sent to the ESP32

Command Sent After Hitting Voice Button


15. Check the Unity terminal to view the results of the voice command. There should only be one command in the console if done correctly. (NOTE: If the Unity engine is on and the phone connects afterwards, there is trash data that is sent from the phone to Unity. To resolve this, simply unplay the Unity engine and hit play again)


Unity console after command received


16. To use the Human Model, flash the ESP32 with the bluetoothApp.ino file and load Unity with the Prosthetic scene found in Assets > Scenes. 


Human Model


17. This model utilizes the use of the ESP32s ADC1 on pins 32, 33, and 34 which correspond to the x, y, and switch axis respectively.


18. Again, use the Android application to connect to the newly flashed ESP32 and hit play on Unity. The voice inputs for this model are only:

Activate Shoulder
Activate Forearm
Activate Wrist
Activate Fingers


19. For this example, we will use the command “Activate Shoulder”. Similar to before, press the voice button and say “Activate Shoulder”. The app will respond by saying “Activating shoulder” and by using the joystick, we are able to move the shoulder component in any direction we want. Continue this method for the next components that need to be moved. (NOTE: The Unity window must be the active window or else the inputs will lag and not show on the screen)


Shoulder of Human Model moved


20. To reset a component, use the space bar on the keyboard to reset the component to its original position


Shoulder reset

