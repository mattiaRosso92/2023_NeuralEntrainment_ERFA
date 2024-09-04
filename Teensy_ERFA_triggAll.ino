/*
#include <Encoder.h>

#include <TouchSwitch.h>
#include <TouchVelocity.h>
#include <TouchVariable.h>
#include <Flicker.h>

#include <Bounce.h>

#include <MIDItouch.h>
#include <MIDIenc.h>
#include <MIDIcontroller.h>
#include <MIDIpot.h>
#include <MIDIdrum.h>
#include <MIDIbutton.h>
*/


/*
 * Pad for finger-tapping w/ trigger to EEG
 * 
 * Portable version with complete Ableton-Teensy interface.
 * On pin18, bnc port to send TTL to EEG.
 * The logging of metronomes' timestamps is enabled.
 * 
 * NB! red pad ranges at rest is ~ 16-24 vs ~180-190 of the blue. Drastic change in threshold multiplier to compensate
 * 
 * The script is pruned to the essential. 
 * 
*/


//MIDI interface (TODO: perhaps remove these, if not used in the script)
#define noteOnCMD               0x90  //MIDI command 144 in integers
#define noteOffCMD              0x80
#define maxVelocity             127
// Define all Events
#define padChannel              1 // tapping logs
#define steadyBeat              2 // standard events
#define phasePlus               3 // perturbations
#define phaseMins               4
#define tempoPlus               5
#define tempoMins               6
// The following are not going into the events
#define tempoPlusStop           7    
#define tempoMinsStop           8  
// This triggers the beginning and the end of EEG epoch (whole trial)
#define eegChannel              13 // Not used in this version

//MIDI note defines for each trigger 
#define padNote              0x00
//MIDI baud rate
#define serialRate           31250
//padS define
#define padPin               15
//BNC
#define bncPin               18
int bncFlag = 0;

//Variables initialize

//for printing
int val=0;
int id=0; //for pads and metronome labels
int thresh=0;
int abovethresh=0; 
int velocity=0;
int go=0;           //flag for start&stop

//for debouncing without delay()
int timeThresh=0;
int lastTime=0;
int previousTime=0;
int currentTime=0;

//Dynamic threshold settings
int n=200;
int sum=0;
int sumerr=0;
double avg;
double sd;


//for amplitude
int valLast=0;   
int onset=0;
int close2peak=0;
//this depends on sensors' sensitivity
int maxPressure = 250;


void setup() 
{
  Serial.begin(serialRate);       //setup serial port for terminal (cause we want to print it to the screen; the value is the speed we communicate with the computer, 9600 is a common value). 31250 is the mIDI serial rate
  pinMode(padPin, INPUT);
  pinMode(bncPin, OUTPUT);
  setThreshold();
  timeThresh=350;           //for debouncing. it just stops one single pad from starting the if loop during x ms (sort of delay() without stopping the program); higher now, also to avoid extra tap as false positives
}



void loop() 
{

  //listen and handle MIDI
  usbMIDI.read();                                                     
  usbMIDI.setHandleNoteOn(djogTrigger);
  
  //store updated time 
  currentTime = millis();                
  
  //waits start trigger to print taps;
  if (go == 0)   //NB: if set to 0, there is no need of any external Start trigger 
  {     
//test
    
    //play pads
    playPad();
  }  

  
}






void djogTrigger (byte channel, byte pitch, byte velocity)    
{

        
  //log metronome timestamps
  if (channel != padChannel  && channel != eegChannel)
  {
    id=channel;
    Serial.print(currentTime);       
    Serial.print(",");          
    Serial.print(id);
    Serial.print(",");
    Serial.println(0);

    if (channel != steadyBeat && bncFlag == 0) 
    {
     digitalWrite(bncPin, HIGH);       
     bncFlag = 1;      
    }
    else if (channel == steadyBeat && bncFlag == 1)
    {
     digitalWrite(bncPin, LOW);       
     bncFlag = 0; 
    }
    
  }

//THIS ONE GOES AWAY
/*
  //control EEG
  if (channel == eegChannel)
  {  
    if (bncFlag == 0)                                                                             
    {
      //single pulse for triggering start in EEG recording; HIGH until the end                
      digitalWrite(bncPin, HIGH);       
      bncFlag = 1;      
      //Serial.print(currentTime);
      //Serial.print(",");
      //Serial.println("Start eeg");//DEBUG
    } 
    else
    {  
      //single pulse for marking END in EEG recording; set to LOW                     
      digitalWrite(bncPin, LOW);       
      bncFlag = 0;      
      //Serial.print(currentTime);
      //Serial.print(",");
      //Serial.println("Stop eeg");//DEBUG
    }   
  }
*/
  
}



void playPad()
{
  
  //read value from pin         
  val = analogRead(padPin); 
 
  //detect onset
  if (val>thresh && abovethresh==0  && currentTime-lastTime>timeThresh)                  
  {                
    id=1;
    onset=currentTime;     //store onset time, to further print once peak is reached
    
    //to Ableton; I had issues sending the note after printing, for conveying velocity. There were missing notes; for the moment, send on onset for demonstrative purposes and eventual data recovery
    midiNoteOn(padNote , maxVelocity , padChannel);   //not full envelope shape;    ****   map velocity onto amplitude

    abovethresh=1;
    
    //update time 
    lastTime=currentTime; 
  }
  
  //detect peak
  if (abovethresh == 1 && val < valLast)  //past the onset, first flex point
  {
    //print      
    Serial.print(onset);       //previously stored onset
    Serial.print(",");         
    Serial.print(id);           //subject ID
    Serial.print(",");
    Serial.println(valLast);       //amplitude (peak)
    //turn OFF; not expressive, yet
    midiNoteOff(padNote , maxVelocity, padChannel);
    //enable next onset detection
    abovethresh=0; 
    
  }                                     
  
  //update amplitude
  valLast = val;

  //map velocity
  velocity = map(valLast, thresh, maxPressure, 1, maxVelocity);
  if (velocity > maxVelocity)
  {
    velocity = maxVelocity;
  }
  
  //turn OFF; here, when it will be expressive
  //midiNoteOff(padNote , maxVelocity, padChannel);
        
}



void setThreshold() 
{                                
  for (int i=0; i<n; i++) 
    {
      sum = sum+analogRead(padPin);
      delay(5);                    
    }
  avg=sum/n;
  
  for (int i=0; i<n; i++) 
    {
      sumerr = sumerr+pow((analogRead(padPin)-avg),2);
      delay(5);                   
    }
  sd=sq((sumerr)/n);
  
  thresh =  avg * 4;  //(avg + 3*sd);    // red pad works best not considering standard deviation; set avg*3 (2/12/20)
} 


void midiNoteOn(byte note, byte midiVelocity, byte channel)
{
  //if(velocity > maxVelocity)
    //velocity = maxVelocity;
  if(velocity > maxVelocity)
    velocity = maxVelocity;
  usbMIDI.sendNoteOn(note, midiVelocity, channel, 0);   //(note, velocity, channel, cable)
}

void midiNoteOff(byte note, byte midiVelocity, byte channel)
{
  usbMIDI.sendNoteOff(note, midiVelocity, channel, 0);
}
