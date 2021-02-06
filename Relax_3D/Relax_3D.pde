import peasy.*;
import java.util.Scanner;
PeasyCam cam;
PFont font;
ArrayList<Bead> beads = new ArrayList<Bead>();
int randomPoints = 50;
boolean action = false;
boolean box = true;
boolean bezOrLines = false;
boolean help = false;
boolean acceleration = false;
boolean showUmd = false;

//Bezier Tube constants
float res = 10.0;

//simulation constants
float alpha = 5; //electrical exponent
float beta = 1; //mechanical exponent
float K = pow(10, 14);//electrical constant
float H = pow(10, -2);//mechanical constant
float dmax = 10;
float dclose = 15;
float Umd = 0.0;

void setup() {
  size(640, 640, P3D);
  surface.setTitle("Knot Relaxation"); 
  cam = new PeasyCam(this, 300);
  cam.setMinimumDistance(1.0);
  cam.setMaximumDistance(700.0);
  font = createFont("Arial", 24);
  reset("random");
}

void draw() {
  //noLoop();
  background(0);
  noFill();
  stroke(255);
  strokeWeight(3);
  lights();
  
  if (action) {
    //int dt = millis();
    //while (millis() - dt < 3) { //maximum is 16.6
    //  update();
    //}
    update();
  }
  if (bezOrLines) {
    makeBezTube();
  } else {
    lines();
  }

  if (box) showBox();
  cam.beginHUD();
  textFont(font);
  //K = pow(10,map(mouseX,0,width,8,15)); 
  //H = pow(10,map(mouseY,height,0,-5,0));
  //text("H=10^"+log(H)/log(10), 20, 55);
  //text("K=10^"+log(K)/log(10), 20, 75);
  if (showUmd) text("U_MD = "+Umd, 10, 30);
  if (help) {
    fill(0); 
    rect(0, 0, width, height);  
    fill(255);
    text("controls:\ndrag mouse to rotate\nscroll to zoom\ndouble click to reset view\nenter : play/pause animation\nb : show/hide box\nr : randomizes beads\nt : toggles bead/bezier view\nh : show/hide this help", 10, 30);
  }
  cam.endHUD();
}

void showBox() {
  stroke(100);
  strokeWeight(2);
  box(200);
}

void lines() {
  if (beads.size() >= 3) {
    for (int i = 0; i < beads.size()-1; i++) { 
      stroke(100);
      strokeWeight(2);
      line(beads.get(i).pos.x, beads.get(i).pos.y, beads.get(i).pos.z, beads.get(i+1).pos.x, beads.get(i+1).pos.y, beads.get(i+1).pos.z);
    }
    for (int i = 0; i <= beads.size()-1; i++) {
      stroke(255);
      strokeWeight(10);
      point(beads.get(i).pos.x, beads.get(i).pos.y, beads.get(i).pos.z);
    }
    stroke(100);
    strokeWeight(2);
    line(beads.get(0).pos.x, beads.get(0).pos.y, beads.get(0).pos.z, beads.get(beads.size()-1).pos.x, beads.get(beads.size()-1).pos.y, beads.get(beads.size()-1).pos.z);
  }
}


void update() {
  Umd = 0;

  for (int i = 0; i < beads.size(); i++) {
    beads.get(i).calculateMechanicalForce(i);
    beads.get(i).calculateElectricalForce(i);
  }

  for (int i = 0; i < beads.size(); i++) {
    beads.get(i).move(i);
  }
}

void makeBezTube() {
  if (beads.size() >= 3) {
    for (int beadIndex = 0; beadIndex < beads.size(); beadIndex++) {
      int nextBead = beadIndex+1;
      int prevBead = beadIndex-1;
      if (beadIndex == beads.size()-1) nextBead = 0;
      if (beadIndex == 0) prevBead = beads.size()-1;
      PVector b0 = beads.get(prevBead).pos;
      PVector b1 = beads.get(beadIndex).pos;
      PVector b2 = beads.get(nextBead).pos;
      PVector prevMidpoint = new PVector(b0.x, b0.y, b0.z);
      prevMidpoint.add(b1);
      prevMidpoint.mult(0.5);
      PVector nextMidpoint = new PVector(b1.x, b1.y, b1.z);
      nextMidpoint.add(b2);
      nextMidpoint.mult(0.5);
      PVector[] p = new PVector[] {
        prevMidpoint, 
        b1, 
        nextMidpoint
      };
      for (float t = 0; t < 1; t += (1/res)) {
        line(
          qP(p[0].x, p[1].x, p[2].x, t), 
          qP(p[0].y, p[1].y, p[2].y, t), 
          qP(p[0].z, p[1].z, p[2].z, t), 
          qP(p[0].x, p[1].x, p[2].x, t+(1/res)), 
          qP(p[0].y, p[1].y, p[2].y, t+(1/res)), 
          qP(p[0].z, p[1].z, p[2].z, t+(1/res))
          );
      }
    }
  }
}

// parametric quadratic bezier interpolation with 
// endpoints p0 and p2 with control point p1
float qP(float p0, float p1, float p2, float t) {
  return p0*(1-2*t+t*t)+2*t*p1-2*t*t*p1+t*t*p2;
}

void reset(String type) {
  if (beads.size() != 0) {
    for (int i = beads.size()-1; i >= 0; i--) {
      beads.remove(i);
    }
  }
  switch(type) {
  case "file":
    try {
      File file = new File("polyknot.txt");
      Scanner input = new Scanner(file);
      int count = 1;
      double x = 0, y = 0, z = 0;
      while (input.hasNextDouble()) {
        switch (count) {
        case 1:
          x = input.nextDouble();
          count++;
          break;
        case 2:
          y = input.nextDouble();
          count++;
          break;
        case 3:
          z = input.nextDouble();
          count = 1;
          beads.add(new Bead((float)(x*10), (float)(y*10), (float)(z*10)));
          break;
        }
      }
      input.close();
    } 
    catch(Exception ex) {
      ex.printStackTrace();
    }
    break;
  case "random":
    //uniform distribution
    for (int i = 0; i < randomPoints; i++) {
      beads.add(new Bead(random(-100, 100), random(-100, 100), random(-100, 100)));
    }
    break;
  case "random lerp":
    beads.add(new Bead(random(-100, 100), random(-100, 100), random(-100, 100)));
    for (int i = 0; i < randomPoints; i+=3) {
      PVector a = new PVector(beads.get(i).pos.x, beads.get(i).pos.y, beads.get(i).pos.z);
      PVector b = new PVector(random(-100, 100), random(-100, 100), random(-100, 100));
      beads.add(new Bead(lerp(a.x, b.x, .33), lerp(a.y, b.y, .33), lerp(a.z, b.z, .33)));
      beads.add(new Bead(lerp(a.x, b.x, .67), lerp(a.y, b.y, .67), lerp(a.z, b.z, .67)));
      beads.add(new Bead(b.x, b.y, b.z));
    }
    break;
  }
}

void keyPressed() {
  if (keyCode == 10 && beads.size() >= 3) { 
    action = !action;
  }
  if (key == 'f') {
    reset("file");
  }
  if (key == 'b') {
    box = !box;
  }
  if (key == 'r') {
    action = false;
    reset("random");
  }
  if (key == 't') {
    bezOrLines = !bezOrLines;
  }
  if (key == 'l') {
    action = false;
    reset("random lerp");
  }
  if (key == 'u') {
    showUmd = !showUmd;
  }
  if (key == 'h') {
    help = !help;
  }
  if (key == '1') {
    acceleration = true;
    alpha = 5; //electrical exponent
    beta = 1; //mechanical exponent
    K = pow(10, -7);//electrical constant
    H = pow(10, -7);//mechanical constant
    dmax = 100;
  }
  if (key == '2') {
    alpha = 5; //electrical exponent
    beta = 1; //mechanical exponent
    K = pow(10, 11);//electrical constant
    H = pow(10, -2);//mechanical constant
    dmax = 0.5;
    dclose = 0.7;
  }
}
