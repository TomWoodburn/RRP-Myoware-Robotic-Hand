const int numReadings = 20;
const int outpin = 7;

int readings1[numReadings];      // the readings from the analog input
int readings2[numReadings];
int readIndex = 0;              // the index of the current reading
int total1 = 0;                  // the running total
int average1 = 0;                // the average

int inputPin1 = A0;

int total2 = 0;                  // the running total
int average2 = 0;                // the average

int inputPin2 = A1;

int threshold1 = 90;
int threshold2 = 90;

void setup() {
  // initialize serial communication with computer:
  Serial.begin(9600);
  pinMode(outpin,OUTPUT);
  // initialize all the readings to 0:
  for (int thisReading = 0; thisReading < numReadings; thisReading++) {
    readings1[thisReading] = 0;
  }
  for (int thisReading = 0; thisReading < numReadings; thisReading++) {
    readings2[thisReading] = 0;
  }
}

void loop() {
  // subtract the last reading:
  total1 = total1 - readings1[readIndex];
  total2 = total2 - readings2[readIndex];
  // read from the sensor:
  readings1[readIndex] = analogRead(inputPin1);
  readings2[readIndex] = analogRead(inputPin2);
  // add the reading to the total:
  total1 = total1 + readings1[readIndex];
  total2 = total2 + readings2[readIndex];
  // advance to the next position in the array:
  readIndex = readIndex + 1;

  // if we're at the end of the array...
  if (readIndex >= numReadings) {
    // ...wrap around to the beginning:
    readIndex = 0;
  }
  
  // calculate the average:
  average1 = total1 / numReadings;
  average2 = total2 / numReadings;

  Serial.print(average1);
  Serial.print(" ");
  Serial.println(average2);
  if (average1 > threshold1 && average2 > threshold2){
    digitalWrite(outpin,HIGH);
  }
  else{
    digitalWrite(outpin,LOW);
  }
  // send it to the computer as ASCII digits
  /*
  Serial.print(average1);
  Serial.print(" ");
  Serial.println(average2);
  */
  delay(1);        // delay in between reads for stability
}
