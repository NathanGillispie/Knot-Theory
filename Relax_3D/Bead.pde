class Bead {
  PVector pos;
  PVector M = new PVector(0, 0, 0);
  PVector E = new PVector(0, 0, 0);
  boolean moveQ = true;

  Bead() {
    pos = new PVector(0, 0, 0);
  }

  Bead(float x_, float y_, float z_) {
    pos = new PVector(x_, y_, z_);
  }

  void calculateMechanicalForce(int beadIndex) {
    if (beadIndex >= beads.size() - 1) beadIndex = -1;
    PVector nextBead = beads.get(beadIndex+1).pos.copy();
    nextBead.sub(this.pos);
    nextBead.setMag(H * pow(nextBead.mag(), beta+1));
    M.add(nextBead);
    beads.get(beadIndex+1).M.sub(nextBead);
  }

  void calculateElectricalForce(int beadIndex) {
    int lastBead = beads.size()-1;
    if (beadIndex == 0) lastBead -= 1;
    for (int i = 2; i <= lastBead - beadIndex; i++) {
      PVector nextBead = beads.get(beadIndex+i).pos.copy();
      nextBead.sub(this.pos);
      nextBead.setMag(-K * pow(nextBead.mag(), -2-alpha));
      E.add(nextBead);
      beads.get(beadIndex+i).E.sub(nextBead);

      if (showUmd) {
        PVector newnextBead = beads.get(beadIndex+1).pos.copy();
        PVector beadI = beads.get(beadIndex+i).pos.copy();
        PVector nextIBead = (beadIndex + i >= beads.size()-1 ? beads.get(0).pos.copy() : beads.get(beadIndex+i+1).pos.copy());
        float MD = pow(distance3D(this.pos, newnextBead, beadI, nextIBead),2);
        newnextBead.sub(this.pos);
        beadI.sub(nextIBead);
        Umd += (newnextBead.mag()*beadI.mag()/MD);
      }
    }
  }

  void move(int beadIndex) {
    if (moveQ) {
      int prevBead = beadIndex - 1;
      int nextBead = beadIndex + 1;
      if (beadIndex == 0)  prevBead = beads.size() - 1;
      if (beadIndex == beads.size() - 1) nextBead = 0;

      for (int i = 0; i < beads.size(); i++) {
        int nextIBead = i + 1;
        if (i >= beads.size() - 1) nextIBead = 0;
        if (i != prevBead && i != beadIndex && i != nextBead) {
          float dist = distance3D(this.pos, beads.get(nextBead).pos, 
            beads.get(i).pos, beads.get(nextIBead).pos);
          if (dist < dclose) {
            moveQ = false;
            beads.get(nextBead).moveQ = false;
            break;
          }
        }
      }

      PVector sum = M.add(E);
      sum.limit(dmax);
      pos.add(sum);
      if (!acceleration) {
        M.set(0, 0);
        E.set(0, 0);
      }
    }

    moveQ = true;
  }

  //code by David Eberly, Geometric Tools https://www.geometrictools.com/
  //licenced under Creative Commons Attribution 4.0 International Licence
  float distance3D(PVector P0, PVector P1, PVector Q0, PVector Q1) {
    float a, b, c, d, e, D, sc, sN, sD, tc, tN, tD;
    float SMALL_NUM = 0.0001;
    PVector u = new PVector(P1.x-P0.x, P1.y-P0.y, P1.z-P0.z);
    PVector v = new PVector(Q1.x-Q0.x, Q1.y-Q0.y, Q1.z-Q0.z);
    PVector w = new PVector(P0.x-Q0.x, P0.y-Q0.y, P0.z-Q0.z);
    a = u.dot(u);       
    b = u.dot(v);
    c = v.dot(v);        
    d = u.dot(w);
    e = v.dot(w);
    D = a*c - b*b;        
    sc = D;
    sN = D;
    sD = D;
    tc = D;
    tN = D;
    tD = D;


    if (D < SMALL_NUM) { 
      sN = 0.0;        
      sD = 1.0;        
      tN = e;
      tD = c;
    } else {               
      sN = (b*e - c*d);
      tN = (a*e - b*d);
      if (sN < 0.0) {        
        sN = 0.0;
        tN = e;
        tD = c;
      } else if (sN > sD) {  
        sN = sD;
        tN = e + b;
        tD = c;
      }
    }

    if (tN < 0.0) {           
      tN = 0.0;
      if (-d < 0.0)
        sN = 0.0;
      else if (-d > a)
        sN = sD;
      else {
        sN = -d;
        sD = a;
      }
    } else if (tN > tD) {     
      tN = tD;
      if ((-d + b) < 0.0)
        sN = 0;
      else if ((-d + b) > a)
        sN = sD;
      else {
        sN = (-d + b);
        sD = a;
      }
    }
    sc = (abs(sN) < SMALL_NUM ? 0.0 : sN / sD);
    tc = (abs(tN) < SMALL_NUM ? 0.0 : tN / tD);

    //float s = sc;
    //float t = tc;

    PVector[] closest = new PVector[2];
    PVector sctP1 = P1.copy();
    sctP1 = sctP1.mult(sc);
    PVector tctQ1 = Q1.copy();
    tctQ1 = tctQ1.mult(tc);
    PVector sctP0 = P0.copy();
    sctP0 = sctP0.mult(1.0-sc);
    PVector tctQ0 = Q0.copy();
    tctQ0 = tctQ0.mult(1.0-tc);
    closest[0] = sctP0.add(sctP1);
    closest[1] = tctQ0.add(tctQ1);
    PVector diff = closest[0].sub(closest[1]);
    return sqrt(diff.dot(diff));
  }
}
