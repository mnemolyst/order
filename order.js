var drawText = false;
var renderInterval = 33; //milliseconds
var gravity = false;
var gravConst = 32.2 * 12 / 1000000 / 8.25; // ft/s^2 -> heads/ms^2 (1 head = 8.25 in)
var framerate = 0;
var hoist = false;
var pause = false;

function crossProd(x, y) {
  return [x[1]*y[2]-x[2]*y[1], x[2]*y[0]-x[0]*y[2], x[0]*y[1]-x[1]*y[0]];
}

function Scene2D(canvas) {
  this.matrix = new Matrix2D();
  this.context = canvas.getContext('2d');
  this.sceneWidth = canvas.width;
  this.sceneHeight = canvas.height;
  this.panX = 0;
  this.panY = 0;
  this.scale = 1.0;
  this.rotation = 0.0;
  this.points = [];
  this.pointsTransformed = [];
  this.pointsLast = [];
  this.pointsAccel = [];
  this.pointsInvMass = [];
  this.numPoints = 0;
  this.items = [];
  this.numItems = 0;
  this.constraints = [];
  this.numConstraints = 0;
  this.speed = 0.1;
  this.constrainOrder = [];

  this.addPoint = function(x, y) {
    this.points[this.points.length] = this.pointsLast[this.pointsLast.length] = x;
    this.points[this.points.length] = this.pointsLast[this.pointsLast.length] = y;

    this.pointsInvMass[this.pointsInvMass.length] = 1;

    return this.numPoints++;
  }

  this.getPoint = function(idx) {
    var idx2D = idx*2;
    return [this.points[idx2D], this.points[idx2D+1]];
  }

  this.addItem = function(item) {
    this.items[this.items.length] = item;

    return this.numItems++;
  }

  this.randomizeConstrainOrder = function() {
    for (var i=0; i<this.items.length; i++) {
      this.constrainOrder[i] = i;
    }

    for (var i=this.constrainOrder.length-1; i>0; i--) {
      var j = Math.floor(Math.random() * (i+1));
      var temp = this.constrainOrder[i];
      this.constrainOrder[i] = this.constrainOrder[j];
      this.constrainOrder[j] = temp;
    }
  }

  this.render = function() {
    var halfWidth = this.sceneWidth*0.5;
    var halfHeight = this.sceneHeight*0.5;

    this.matrix.identity();
    this.matrix.translate(halfWidth-this.panX, halfHeight-this.panY);
    this.matrix.scale(this.scale, this.scale);
    //this.matrix.rotate(this.rotation);

    this.pointsTransformed = this.matrix.transformArray(this.points);

    this.context.save();
    this.context.fillStyle = '#000000';
    this.context.fillRect(0, 0, this.sceneWidth, this.sceneHeight);
    if (drawText) {
      this.context.fillStyle = '#aaaaaa';
      this.context.textBaseline = 'top';
      this.context.fillText('rotation: ' + this.rotation*180/Math.PI, 5, 5);
      this.context.fillText('framerate: ' + framerate, 5, 20);
    }

    for (var i=0; i<this.numItems; i++) {
      this.items[i].render(this);
    }

    this.context.restore();
  }

  this.accumulateForces = function() {
    for (var i=0; i<this.numPoints; i++) {
      var i2 = i*2;
      var xIdx = i2;
      var yIdx = i2+1;

      //this.points3DAccel[xIdx] = gravity?gravConst*windCos*0.1:0;
      //this.points3DAccel[yIdx] = gravity?gravConst*windSin*0.1:0;
      this.pointsAccel[xIdx] = 0;
      this.pointsAccel[yIdx] = gravity?-gravConst:0; // 1g of gravity in heads (8.25in: female) per millisec^2
    }

    //TODO poll each soft constraint item to get its accel contribs
  }

  this.verlet = function(timeDel) {
    for (var i=0; i<this.numPoints; i++) {
      var i2 = i*2;
      var xIdx = i2;
      var yIdx = i2+1;

      var tempX = this.points[xIdx];
      var tempY = this.points[yIdx];

      this.points[xIdx] += this.points[xIdx] - this.pointsLast[xIdx] + this.pointsAccel[xIdx]*timeDel*timeDel;
      this.points[yIdx] += this.points[yIdx] - this.pointsLast[yIdx] + this.pointsAccel[yIdx]*timeDel*timeDel;

      this.pointsLast[xIdx] = tempX;
      this.pointsLast[yIdx] = tempY;
    }
  }

  this.constrain = function() {
    for (var num=0; num<1; num++) {
      for (var i=0; i<this.numItems; i++) {
        this.items[i].constrain();
      }
    }
  }

  this.clip = function() {
    var temp = 0;
    for (var i=0; i<this.numPoints; i++) {
      if (this.points[i*2] < -this.sceneWidth*0.5 || this.points[i*2] > this.sceneWidth*0.5) {
        temp = this.points[i*2];
        this.points[i*2] = this.pointsLast[i*2];
        this.pointsLast[i*2] = temp;
      }
      if (this.points[i*2+1] < -this.sceneHeight*0.5 || this.points[i*2+1] > this.sceneHeight*0.5) {
        temp = this.points[i*2+1];
        this.points[i*2+1] = this.pointsLast[i*2+1];
        this.pointsLast[i*2+1] = temp;
      }
    }
  }

  this.collide = function() {
    for (var i=0; i<this.numItems; i++) {
      this.items[i].color = '#eeeeee';
    }
    for (var i=0; i<this.numItems; i++) {
      for (var j=i+1; j<this.numItems; j++) {
        if (this.items[i].collide(this.items[j])) {
          this.items[i].color = '#ee0000';
          this.items[j].color = '#ee0000';
        }
      }
    }
  }

  this.tic = function() {
    this.render();
    this.accumulateForces();
    this.verlet(renderInterval * this.speed);
    this.constrain();
    //this.collide();
    //this.collide();
    this.constrain();
    //this.collide();
    //this.collide();
    this.clip();
  }
}

function polyItem(scene, pIdx) {
  this.scene = scene;
  this.pIdx = pIdx;
  this.p2Didx = [];
  for (var i=0; i<this.pIdx.length; i++) {
    this.p2Didx[i] = 2 * this.pIdx[i];
  }
  this.sideLengths = [];
  this.restLengths = [];
  this.restLengthsSq = [];
  this.edgeVecIdx = [];
  this.normalDepths = [];
  this.color = '#eeeeee';

  // INIT
  var points = this.scene.points;
  var minSelfProj = 0;
  for (var i=0; i<this.p2Didx.length; i++) {
    var si = (i+1<this.p2Didx.length)?i+1:0; // successor of i

    // compute side lengths
    var v = [points[this.p2Didx[si]]-points[this.p2Didx[i]], points[this.p2Didx[si]+1]-points[this.p2Didx[i]+1]];
    this.sideLengths[i] = Math.sqrt(v[0]*v[0] + v[1]*v[1]);

    // constraint length pairs
    this.edgeVecIdx[this.edgeVecIdx.length] = this.p2Didx[i];
    this.edgeVecIdx[this.edgeVecIdx.length] = this.p2Didx[si];
    var n = this.p2Didx.length-3-(i>1?i-1:0);
    for (var j=i+2; j<i+n+2; j++) {
      this.edgeVecIdx[this.edgeVecIdx.length] = this.p2Didx[i];
      this.edgeVecIdx[this.edgeVecIdx.length] = this.p2Didx[j];
    }

    // depth from each face
    var norml = [(points[this.p2Didx[si]+1]-points[this.p2Didx[i]+1])/this.sideLengths[i], // right-hand unit normal of i-th side
                 (points[this.p2Didx[i]]-points[this.p2Didx[si]])/this.sideLengths[i]];
    this.normalDepths[i] = 0;
    for (var j=0; j<this.p2Didx.length; j++) {
      var proj = (points[this.p2Didx[j]]-points[this.p2Didx[i]])*norml[0] + (points[this.p2Didx[j]+1]-points[this.p2Didx[i]+1])*norml[1];
      this.normalDepths[i] = (this.normalDepths[i] < -proj)?-proj:this.normalDepths[i];
    }
  }

  // constraint rest lengths
  for (var i=0; i<this.edgeVecIdx.length/2; i++) {
    var v = [points[this.edgeVecIdx[i*2+1]]-points[this.edgeVecIdx[i*2]], points[this.edgeVecIdx[i*2+1]+1]-points[this.edgeVecIdx[i*2]+1]];
    this.restLengths[i] = Math.sqrt(v[0]*v[0] + v[1]*v[1]);
    this.restLengthsSq[i] = this.restLengths[i]*this.restLengths[i];
  }

  // METHODS
  this.render = function() {
    var points = this.scene.pointsTransformed;
    var context = this.scene.context;

    context.beginPath();
    context.strokeStyle = this.color;
    for (var i=0; i<this.p2Didx.length; i++) {
      context.lineTo(points[this.p2Didx[i]], this.scene.sceneHeight-points[this.p2Didx[i]+1]);
    }
    context.lineTo(points[this.p2Didx[0]], this.scene.sceneHeight-points[this.p2Didx[0]+1]); // close the polygon
    context.stroke();
  }

  this.constrain = function() {
    var points = this.scene.points;
    for (var i=0; i<this.edgeVecIdx.length/2; i++) { //TODO 
      var vecH = i*2;
      var delta = [points[this.edgeVecIdx[vecH+1]]-points[this.edgeVecIdx[vecH]],
                   points[this.edgeVecIdx[vecH+1]+1]-points[this.edgeVecIdx[vecH]+1]];
      var delSq = delta[0]*delta[0] + delta[1]*delta[1];
      var im1 = 1;//this.scene.pointsInvMass[this.edgeVecIdx[vecH]/2];
      var im2 = 1;//this.scene.pointsInvMass[this.edgeVecIdx[vecH+1]/2];
      var diff = (delSq-this.restLengthsSq[i])/((delSq+this.restLengthsSq[i])*(im1+im2)*2);
      for (var j=0; j<2; j++) {
        var del = diff*delta[j];
        points[this.edgeVecIdx[vecH]+j]   += del*im1;
        points[this.edgeVecIdx[vecH+1]+j] -= del*im2;
      }
    }
  }

  this.collide = function(poly2) {
    var points = this.scene.points;
    //var thisCenter = [points[this.p2Didx[0]]+points[this.p2Didx[1]]+points[this.p2Didx[2]]+points[this.p2Didx[3]],
    //                  points[this.p2Didx[0]+1]+points[this.p2Didx[1]+1]+points[this.p2Didx[2]+1]+points[this.p2Didx[3]+1]];
    //var thatCenter = [(points[poly2.p2Didx[0]]+points[poly2.p2Didx[1]]+points[poly2.p2Didx[2]]+points[poly2.p2Didx[3]])/4.0,
    //                  (points[poly2.p2Didx[0]+1]+points[poly2.p2Didx[1]+1]+points[poly2.p2Didx[2]+1]+points[poly2.p2Didx[3]+1])/4.0];

    var min12Overlap = 0;
    var min21Overlap = 0;
    var min12Norml = [];
    var min21Norml = [];
    var besti = -1;
    var bestj = -1;
    var bestp = -1;
    var doesOverlap = true;

    // First check this poly's sides for penetration against "that" poly's points
    for (var i=0; i<this.p2Didx.length; i++) {
      var j = (i+1<this.p2Didx.length)?i+1:0;
      var norml = [(points[this.p2Didx[j]+1]-points[this.p2Didx[i]+1])/this.restLengths[i], // unit normal of i-th side
                   (points[this.p2Didx[i]]-points[this.p2Didx[j]])/this.restLengths[i]];

      var min12Proj = 0.0;
      var max12Proj = 0.0;
      for (var p=0; p<poly2.p2Didx.length; p++) {
        var proj = (points[poly2.p2Didx[p]]-points[this.p2Didx[i]])*norml[0] + (points[poly2.p2Didx[p]+1]-points[this.p2Didx[i]+1])*norml[1];
        if (proj<min12Proj) {
          min12Proj = proj;
          min12Projp = p;
        }
        max12Proj = (max12Proj<proj)?proj:max12Proj;
      }

      var minProj = (-this.normalDepths[i]<min12Proj)?-this.normalDepths[i]:min12Proj;
      var maxProj = (max12Proj>0.0)?max12Proj:0.0;

      var overlap = max12Proj - min12Proj + this.normalDepths[i] - maxProj + minProj;
      if (overlap > 0) {
        doesOverlap = doesOverlap && true;
        if (overlap < min12Overlap) {
          min12Overlap = overlap;
          besti = i;
          bestj = j;
          bestp = p;
        }
      }

      if (!doesOverlap)
        return false;
    }

    // Then check "that" poly's sides
    for (var i=0; i<4; i++) {
      var j = i>2?i-3:i+1; // wrap around
      var k = i>1?i-2:i+2;
      var l = i>0?i-1:i+3;
      var norml = [(points[poly2.p2Didx[j]+1]-points[poly2.p2Didx[i]+1])/poly2.restLengths[i], // unit normal of i-th side
                   (points[poly2.p2Didx[i]]-points[poly2.p2Didx[j]])/poly2.restLengths[i]];

      var proj2a = (points[poly2.p2Didx[k]]-points[poly2.p2Didx[i]])*norml[0] + (points[poly2.p2Didx[k]+1]-points[poly2.p2Didx[i]+1])*norml[1];
      var proj2b = (points[poly2.p2Didx[l]]-points[poly2.p2Didx[i]])*norml[0] + (points[poly2.p2Didx[l]+1]-points[poly2.p2Didx[i]+1])*norml[1];
      var min2Proj = (proj2a<proj2b)?proj2a:proj2b;

      var min1Proj = 0.0;
      var max1Proj = 0.0;
      for (var m=0; m<4; m++) {
        var proj = (points[this.p2Didx[m]]-points[poly2.p2Didx[i]])*norml[0] + (points[this.p2Didx[m]+1]-points[poly2.p2Didx[i]+1])*norml[1];
        min1Proj = (proj<min1Proj)?proj:min1Proj;
        max1Proj = (max1Proj<proj)?proj:max1Proj;

        if (min1Proj < 0 && min1Proj > maxNeg1Proj) {
          maxNeg1Proj = min1Proj;
          maxNeg1Projm = m;
          maxNeg1ProjNorml = norml;
          maxNeg1Proji = i;
          maxNeg1Projj = j;
        }
      }

      var minProj = (min1Proj<min2Proj)?min1Proj:min2Proj;
      var maxProj = (max1Proj>0.0)?max1Proj:0.0;

      var overlaps = overlaps && ((max1Proj - min1Proj - min2Proj) > (maxProj - minProj));

      if (!overlaps)
        return false;
    }

    // Update points
    if (maxNeg1Proj > maxNeg2Proj) {    // a point on this poly hit a side on that poly
      if (Math.abs(points[poly2.p2Didx[maxNeg1Projj]] - points[poly2.p2Didx[maxNeg1Proji]]) > Math.abs(points[poly2.p2Didx[maxNeg1Projj]+1] - points[poly2.p2Didx[maxNeg1Proji]+1]))
        var ijHitF = Math.abs(points[this.p2Didx[maxNeg1Projm]]-points[poly2.p2Didx[maxNeg1Proji]]) / Math.abs(points[poly2.p2Didx[maxNeg1Projj]]-points[poly2.p2Didx[maxNeg1Proji]]);
      else
        var ijHitF = Math.abs(points[this.p2Didx[maxNeg1Projm]+1]-points[poly2.p2Didx[maxNeg1Proji]+1]) / Math.abs(points[poly2.p2Didx[maxNeg1Projj]+1]-points[poly2.p2Didx[maxNeg1Proji]+1]);
      var l = 1.0/(ijHitF*ijHitF + (1-ijHitF)*(1-ijHitF));
      var iHitP = (1-ijHitF)*l;
      var jHitP = ijHitF*l;
      var push = [maxNeg1ProjNorml[0]*maxNeg1Proj*0.5, maxNeg1ProjNorml[1]*maxNeg1Proj*0.5];
      points[poly2.p2Didx[maxNeg1Proji]]   += push[0]*iHitP;
      points[poly2.p2Didx[maxNeg1Proji]+1] += push[1]*iHitP;
      points[poly2.p2Didx[maxNeg1Projj]]   += push[0]*jHitP;
      points[poly2.p2Didx[maxNeg1Projj]+1] += push[1]*jHitP;
      points[this.p2Didx[maxNeg1Projm]]    -= push[0];
      points[this.p2Didx[maxNeg1Projm]+1]  -= push[1];
    } else {                            // a _side_ on this poly hit a _point_ on that poly
      if (Math.abs(points[this.p2Didx[maxNeg2Projj]] - points[this.p2Didx[maxNeg2Proji]]) > Math.abs(points[this.p2Didx[maxNeg2Projj]+1]-points[this.p2Didx[maxNeg2Proji]+1]))
        var ijHitF = Math.abs(points[poly2.p2Didx[maxNeg2Projm]]-points[this.p2Didx[maxNeg2Proji]]) / Math.abs(points[this.p2Didx[maxNeg2Projj]]-points[this.p2Didx[maxNeg2Proji]]);
      else
        var ijHitF = Math.abs(points[poly2.p2Didx[maxNeg2Projm]+1]-points[this.p2Didx[maxNeg2Proji]+1]) / Math.abs(points[this.p2Didx[maxNeg2Projj]+1]-points[this.p2Didx[maxNeg2Proji]+1]);
      var l = 1.0/(ijHitF*ijHitF + (1-ijHitF)*(1-ijHitF));
      var iHitP = (1-ijHitF)*l;
      var jHitP = ijHitF*l;
      var push = [maxNeg2ProjNorml[0]*maxNeg2Proj*0.5, maxNeg2ProjNorml[1]*maxNeg2Proj*0.5];
      points[this.p2Didx[maxNeg2Proji]]    += push[0]*iHitP;
      points[this.p2Didx[maxNeg2Proji]+1]  += push[1]*iHitP;
      points[this.p2Didx[maxNeg2Projj]]    += push[0]*jHitP;
      points[this.p2Didx[maxNeg2Projj]+1]  += push[1]*jHitP;
      points[poly2.p2Didx[maxNeg2Projm]]   -= push[0];
      points[poly2.p2Didx[maxNeg2Projm]+1] -= push[1];
    }

    return true;
  }
}

function quadItem(scene, p1idx, p2idx, p3idx, p4idx) {
  this.scene = scene;
  this.p1idx = p1idx;
  this.p2idx = p2idx;
  this.p3idx = p3idx;
  this.p4idx = p4idx;
  this.pIdx = [
    p1idx,
    p2idx,
    p3idx,
    p4idx,];
  this.p2Didx = [
    p1idx*2,
    p2idx*2,
    p3idx*2,
    p4idx*2,]; // p1 repeated for circularity
  this.color = '#eeeeee';

  // INIT
  // compute rest lengths etc
  var p = [
    scene.getPoint(p1idx),
    scene.getPoint(p2idx),
    scene.getPoint(p3idx),
    scene.getPoint(p4idx),];
  this.edgeVecIdx = [
    0,1,
    1,2,
    2,3,
    3,0,
    0,2,
    1,3,];
  this.restLengths = [];
  this.restLengthsSq = [];
  for (var i=0; i<6; i++) {
    var rad = (p[this.edgeVecIdx[i*2+1]][0]-p[this.edgeVecIdx[i*2]][0])*(p[this.edgeVecIdx[i*2+1]][0]-p[this.edgeVecIdx[i*2]][0]);
        rad +=(p[this.edgeVecIdx[i*2+1]][1]-p[this.edgeVecIdx[i*2]][1])*(p[this.edgeVecIdx[i*2+1]][1]-p[this.edgeVecIdx[i*2]][1]);
    this.restLengths[i] = Math.sqrt(rad);
    this.restLengthsSq[i] = this.restLengths[i]*this.restLengths[i];
  }
  this.overlappedSides = [false, false, false, false];

  // METHODS
  this.render = function() {
    var points = this.scene.pointsTransformed;
    var context = this.scene.context;
    var path = [0,1,2,3,0];

    context.beginPath();
    context.strokeStyle = this.color;
    for (var i=0; i<5; i++) {
      //alert(points[this.p2Didx[path[i]]] + ' ' + points[this.p2Didx[path[i]]+1]);
      context.lineTo(points[this.p2Didx[path[i]]], this.scene.sceneHeight-points[this.p2Didx[path[i]]+1]);
    }
    context.stroke();
  }

  this.constrain = function() {
    var points = this.scene.points;
    for (var i=0; i<6; i++) {
      var vecH = i*2;
      var delta = [points[this.p2Didx[this.edgeVecIdx[vecH+1]]]-points[this.p2Didx[this.edgeVecIdx[vecH]]],
                   points[this.p2Didx[this.edgeVecIdx[vecH+1]]+1]-points[this.p2Didx[this.edgeVecIdx[vecH]]+1]];
      var delSq = delta[0]*delta[0] + delta[1]*delta[1];
      var im1 = this.scene.pointsInvMass[this.pIdx[this.edgeVecIdx[vecH]]];
      var im2 = this.scene.pointsInvMass[this.pIdx[this.edgeVecIdx[vecH+1]]];
      var diff = (delSq-this.restLengthsSq[i])/((delSq+this.restLengthsSq[i])*(im1+im2)*2);
      for (var j=0; j<2; j++) {
        var del = diff*delta[j];
        points[this.p2Didx[this.edgeVecIdx[vecH]]+j]   += del*im1;
        points[this.p2Didx[this.edgeVecIdx[vecH+1]]+j] -= del*im2;
      }
    }
  }

  this.collide = function(quad2) {
    var points = this.scene.points;
    var thisCenter = [points[this.p2Didx[0]]+points[this.p2Didx[1]]+points[this.p2Didx[2]]+points[this.p2Didx[3]],
                      points[this.p2Didx[0]+1]+points[this.p2Didx[1]+1]+points[this.p2Didx[2]+1]+points[this.p2Didx[3]+1]];
    var thatCenter = [(points[quad2.p2Didx[0]]+points[quad2.p2Didx[1]]+points[quad2.p2Didx[2]]+points[quad2.p2Didx[3]])/4.0,
                      (points[quad2.p2Didx[0]+1]+points[quad2.p2Didx[1]+1]+points[quad2.p2Didx[2]+1]+points[quad2.p2Didx[3]+1])/4.0];

    var maxNeg1Proj = -1000;
    var maxNeg2Proj = -1000;
    var maxNeg1Projm = 0;
    var maxNeg2Projm = 0;
    var maxNeg1ProjNorml = [];
    var maxNeg2ProjNorml = [];
    var maxNeg1Proji = -1;
    var maxNeg1Projj = -1;
    var maxNeg1Projm = -1;
    var maxNeg2Proji = -1;
    var maxNeg2Projj = -1;
    var maxNeg2Projm = -1;
    var overlaps = true;

    // First check this quad's sides for penetration against "that" quad's points
    for (var i=0; i<4; i++) {
      var j = i>2?i-3:i+1; // wrap around
      var k = i>1?i-2:i+2;
      var l = i>0?i-1:i+3;
      var norml = [(points[this.p2Didx[j]+1]-points[this.p2Didx[i]+1])/this.restLengths[i], // unit normal of i-th side
                   (points[this.p2Didx[i]]-points[this.p2Didx[j]])/this.restLengths[i]];

      // which of this quad's points is furthest "left"?
      var proj1a = (points[this.p2Didx[k]]-points[this.p2Didx[i]])*norml[0] + (points[this.p2Didx[k]+1]-points[this.p2Didx[i]+1])*norml[1];
      var proj1b = (points[this.p2Didx[l]]-points[this.p2Didx[i]])*norml[0] + (points[this.p2Didx[l]+1]-points[this.p2Didx[i]+1])*norml[1];
      var min1Proj = (proj1a<proj1b)?proj1a:proj1b;

      var min2Proj = 0.0;
      var max2Proj = 0.0;
      for (var m=0; m<4; m++) {
        var proj = (points[quad2.p2Didx[m]]-points[this.p2Didx[i]])*norml[0] + (points[quad2.p2Didx[m]+1]-points[this.p2Didx[i]+1])*norml[1];
        min2Proj = (proj<min2Proj)?proj:min2Proj;
        max2Proj = (max2Proj<proj)?proj:max2Proj;

        if (min2Proj < 0 && min2Proj > maxNeg2Proj) {
          maxNeg2Proj = min2Proj;
          maxNeg2Projm = m;
          maxNeg2ProjNorml = norml;
          maxNeg2Proji = i;
          maxNeg2Projj = j;
        }
      }

      var minProj = (min1Proj<min2Proj)?min1Proj:min2Proj;
      var maxProj = (max2Proj>0.0)?max2Proj:0.0;

      var overlaps = overlaps && ((max2Proj - min2Proj - min1Proj) > (maxProj - minProj));

      if (!overlaps)
        return false;
    }

    // Then check "that" quad's sides
    for (var i=0; i<4; i++) {
      var j = i>2?i-3:i+1; // wrap around
      var k = i>1?i-2:i+2;
      var l = i>0?i-1:i+3;
      var norml = [(points[quad2.p2Didx[j]+1]-points[quad2.p2Didx[i]+1])/quad2.restLengths[i], // unit normal of i-th side
                   (points[quad2.p2Didx[i]]-points[quad2.p2Didx[j]])/quad2.restLengths[i]];

      var proj2a = (points[quad2.p2Didx[k]]-points[quad2.p2Didx[i]])*norml[0] + (points[quad2.p2Didx[k]+1]-points[quad2.p2Didx[i]+1])*norml[1];
      var proj2b = (points[quad2.p2Didx[l]]-points[quad2.p2Didx[i]])*norml[0] + (points[quad2.p2Didx[l]+1]-points[quad2.p2Didx[i]+1])*norml[1];
      var min2Proj = (proj2a<proj2b)?proj2a:proj2b;

      var min1Proj = 0.0;
      var max1Proj = 0.0;
      for (var m=0; m<4; m++) {
        var proj = (points[this.p2Didx[m]]-points[quad2.p2Didx[i]])*norml[0] + (points[this.p2Didx[m]+1]-points[quad2.p2Didx[i]+1])*norml[1];
        min1Proj = (proj<min1Proj)?proj:min1Proj;
        max1Proj = (max1Proj<proj)?proj:max1Proj;

        if (min1Proj < 0 && min1Proj > maxNeg1Proj) {
          maxNeg1Proj = min1Proj;
          maxNeg1Projm = m;
          maxNeg1ProjNorml = norml;
          maxNeg1Proji = i;
          maxNeg1Projj = j;
        }
      }

      var minProj = (min1Proj<min2Proj)?min1Proj:min2Proj;
      var maxProj = (max1Proj>0.0)?max1Proj:0.0;

      var overlaps = overlaps && ((max1Proj - min1Proj - min2Proj) > (maxProj - minProj));

      if (!overlaps)
        return false;
    }

    // Update points
    if (maxNeg1Proj > maxNeg2Proj) {    // a point on this quad hit a side on that quad
      if (Math.abs(points[quad2.p2Didx[maxNeg1Projj]] - points[quad2.p2Didx[maxNeg1Proji]]) > Math.abs(points[quad2.p2Didx[maxNeg1Projj]+1] - points[quad2.p2Didx[maxNeg1Proji]+1]))
        var ijHitF = Math.abs(points[this.p2Didx[maxNeg1Projm]]-points[quad2.p2Didx[maxNeg1Proji]]) / Math.abs(points[quad2.p2Didx[maxNeg1Projj]]-points[quad2.p2Didx[maxNeg1Proji]]);
      else
        var ijHitF = Math.abs(points[this.p2Didx[maxNeg1Projm]+1]-points[quad2.p2Didx[maxNeg1Proji]+1]) / Math.abs(points[quad2.p2Didx[maxNeg1Projj]+1]-points[quad2.p2Didx[maxNeg1Proji]+1]);
      var l = 1.0/(ijHitF*ijHitF + (1-ijHitF)*(1-ijHitF));
      var iHitP = (1-ijHitF)*l;
      var jHitP = ijHitF*l;
      var push = [maxNeg1ProjNorml[0]*maxNeg1Proj*0.5, maxNeg1ProjNorml[1]*maxNeg1Proj*0.5];
      points[quad2.p2Didx[maxNeg1Proji]]   += push[0]*iHitP;
      points[quad2.p2Didx[maxNeg1Proji]+1] += push[1]*iHitP;
      points[quad2.p2Didx[maxNeg1Projj]]   += push[0]*jHitP;
      points[quad2.p2Didx[maxNeg1Projj]+1] += push[1]*jHitP;
      points[this.p2Didx[maxNeg1Projm]]    -= push[0];
      points[this.p2Didx[maxNeg1Projm]+1]  -= push[1];
    } else {                            // a _side_ on this quad hit a _point_ on that quad
      if (Math.abs(points[this.p2Didx[maxNeg2Projj]] - points[this.p2Didx[maxNeg2Proji]]) > Math.abs(points[this.p2Didx[maxNeg2Projj]+1]-points[this.p2Didx[maxNeg2Proji]+1]))
        var ijHitF = Math.abs(points[quad2.p2Didx[maxNeg2Projm]]-points[this.p2Didx[maxNeg2Proji]]) / Math.abs(points[this.p2Didx[maxNeg2Projj]]-points[this.p2Didx[maxNeg2Proji]]);
      else
        var ijHitF = Math.abs(points[quad2.p2Didx[maxNeg2Projm]+1]-points[this.p2Didx[maxNeg2Proji]+1]) / Math.abs(points[this.p2Didx[maxNeg2Projj]+1]-points[this.p2Didx[maxNeg2Proji]+1]);
      var l = 1.0/(ijHitF*ijHitF + (1-ijHitF)*(1-ijHitF));
      var iHitP = (1-ijHitF)*l;
      var jHitP = ijHitF*l;
      var push = [maxNeg2ProjNorml[0]*maxNeg2Proj*0.5, maxNeg2ProjNorml[1]*maxNeg2Proj*0.5];
      points[this.p2Didx[maxNeg2Proji]]    += push[0]*iHitP;
      points[this.p2Didx[maxNeg2Proji]+1]  += push[1]*iHitP;
      points[this.p2Didx[maxNeg2Projj]]    += push[0]*jHitP;
      points[this.p2Didx[maxNeg2Projj]+1]  += push[1]*jHitP;
      points[quad2.p2Didx[maxNeg2Projm]]   -= push[0];
      points[quad2.p2Didx[maxNeg2Projm]+1] -= push[1];
    }

    return true;
  }
}

function linItem(scene, p1idx, p2idx) {
  this.scene = scene;
  this.p1idx = p1idx;
  this.p2idx = p2idx;
  this.p12Didx = p1idx*2;
  this.p22Didx = p2idx*2;
  this.p13Didx = p1idx*3;
  this.p23Didx = p2idx*3;

  // INIT
  // compute rest length
  var p1 = scene.getPoint3D(p1idx);
  var p2 = scene.getPoint3D(p2idx);
  this.restLength = Math.sqrt((p2[0]-p1[0])*(p2[0]-p1[0]) + (p2[1]-p1[1])*(p2[1]-p1[1]) + (p2[2]-p1[2])*(p2[2]-p1[2]));
  this.restSq = this.restLength * this.restLength;

  // METHODS
  this.render = function() {
    var points = this.scene.points2D;
    var context = this.scene.context;

    context.beginPath();
    context.strokeStyle = '#eeeeee';
    context.lineTo(points[this.p12Didx], points[this.p12Didx+1]);
    context.lineTo(points[this.p22Didx], points[this.p22Didx+1]);
    context.stroke();
  }

  this.constrain= function() {
    var points = this.scene.points3D;
    var delta = [points[this.p23Didx]-points[this.p13Didx], points[this.p23Didx+1]-points[this.p13Didx+1], points[this.p23Didx+2]-points[this.p13Didx+2]];
    var delSq = delta[0]*delta[0] + delta[1]*delta[1] + delta[2]*delta[2];
    var diff = (delSq-this.restSq)/(delSq+this.restSq)*0.5;
    for (var j=0; j<3; j++) {
      points[this.p13Didx+j] += diff * delta[j];
      points[this.p23Didx+j] -= diff * delta[j];
    }
  }
}

function cordItem(scene, maxLength, points) {
  this.scene = scene;
  this.maxLength = maxLength;
  this.avgMaxSeg = maxLength / (points.length - 1);
  this.points = points;
  this.p2Didx = [];
  this.p3Didx = [];
  for (var i=0; i<points.length; i++) {
    this.p2Didx[i] = points[i]*2;
    this.p3Didx[i] = points[i]*3;
  }

  // METHODS
  this.render = function() {
    var points = this.scene.points2D;
    var context = this.scene.context;

    context.beginPath();
    context.strokeStyle = '#ee6666';
    for (var i=0; i<this.p2Didx.length; i++) {
      context.lineTo(points[this.p2Didx[i]], points[this.p2Didx[i]+1]);
    }
    context.stroke();
  }

  this.constrain= function() {
    var points = this.scene.points3D;
    var deltas = [];
    var delSqs = [];
    var delLens = [];
    var totalLen = 0;
    for (var i=1; i<this.p3Didx.length; i++) {
      var delta = deltas[i-1] = [points[this.p3Didx[i]]-points[this.p3Didx[i-1]],
                                 points[this.p3Didx[i]+1]-points[this.p3Didx[i-1]+1],
                                 points[this.p3Didx[i]+2]-points[this.p3Didx[i-1]+2]];
      delSqs[i-1] = delta[0]*delta[0] + delta[1]*delta[1] + delta[2]*delta[2];
      //totalLen += delLens[i-1] = Math.sqrt(delSqs[i-1]);
      var apx = delSqs[i-1];
      apx = (apx + delSqs[i-1]/apx)*0.5;
      apx = (apx + delSqs[i-1]/apx)*0.5;
      apx = (apx + delSqs[i-1]/apx)*0.5;
      totalLen += delLens[i-1] = apx;
    }
    if (totalLen <= maxLength)
      return;

    var err = (totalLen-maxLength) / (this.p3Didx.length-1);
    for (var i=1; i<this.p3Didx.length; i++) {
      var diff = err/(delLens[i-1]*2);
      for (var j=0; j<3; j++) {
        points[this.p3Didx[i-1]+j] += diff * deltas[i-1][j];
        points[this.p3Didx[i]+j] -= diff * deltas[i-1][j];
      }
    }
  }
}

function linSumPairItem(scene, p1idx, p2idx, p3idx, p4idx, maxShrink, maxGrow) {
  this.scene = scene;
  this.p1idx = p1idx;
  this.p2idx = p2idx;
  this.p3idx = p3idx;
  this.p4idx = p4idx;
  this.maxShrink = maxShrink;
  this.maxGrow = maxGrow;
  this.p12Didx = p1idx*2;
  this.p22Didx = p2idx*2;
  this.p32Didx = p3idx*2;
  this.p42Didx = p4idx*2;
  this.p13Didx = p1idx*3;
  this.p23Didx = p2idx*3;
  this.p33Didx = p3idx*3;
  this.p43Didx = p4idx*3;

  // INIT
  // compute rest length
  var p1 = scene.getPoint3D(p1idx);
  var p2 = scene.getPoint3D(p2idx);
  var p3 = scene.getPoint3D(p3idx);
  var p4 = scene.getPoint3D(p4idx);
  this.restLength = Math.sqrt((p2[0]-p1[0])*(p2[0]-p1[0]) + (p2[1]-p1[1])*(p2[1]-p1[1]) + (p2[2]-p1[2])*(p2[2]-p1[2]));
  this.restLength += Math.sqrt((p4[0]-p3[0])*(p4[0]-p3[0]) + (p4[1]-p3[1])*(p4[1]-p3[1]) + (p4[2]-p3[2])*(p4[2]-p3[2]));
  this.restSq = this.restLength * this.restLength;
  this.minLength = this.restLength - maxShrink;
  this.maxLength = this.restLength + maxGrow;

  // METHODS
  this.render = function() {
    var points = this.scene.points2D;
    var context = this.scene.context;

    context.beginPath();
    context.strokeStyle = '#ee3333';
    context.lineTo(points[this.p12Didx], points[this.p12Didx+1]);
    context.lineTo(points[this.p22Didx], points[this.p22Didx+1]);
    context.stroke();
    context.moveTo(points[this.p32Didx], points[this.p32Didx+1]);
    context.lineTo(points[this.p42Didx], points[this.p42Didx+1]);
    context.stroke();
  }

  this.constrain= function() {
    var points = this.scene.points3D;
    var delta1 = [points[this.p23Didx]-points[this.p13Didx], points[this.p23Didx+1]-points[this.p13Didx+1], points[this.p23Didx+2]-points[this.p13Didx+2]];
    var delta2 = [points[this.p43Didx]-points[this.p33Didx], points[this.p43Didx+1]-points[this.p33Didx+1], points[this.p43Didx+2]-points[this.p33Didx+2]];
    var delSq1 = delta1[0]*delta1[0] + delta1[1]*delta1[1] + delta1[2]*delta1[2];
    var delSq2 = delta2[0]*delta2[0] + delta2[1]*delta2[1] + delta2[2]*delta2[2];
    var delLen = Math.sqrt(delSq1) + Math.sqrt(delSq2); //TODO optimize

    if (delLen > this.maxLength) {
      var diff = (delLen - this.maxLength) / delLen * 0.23;
    } else if (delLen < this.minLength) {
      var diff = (delLen - this.minLength) / delLen * 0.23;
    } else {
      return;
    }

    for (var j=0; j<3; j++) {
      points[this.p13Didx+j] += diff * delta1[j];
      points[this.p23Didx+j] -= diff * delta1[j];
      points[this.p33Didx+j] += diff * delta2[j];
      points[this.p43Didx+j] -= diff * delta2[j];
    }
  }
}

function linDiffPairItem(scene, p1idx, p2idx, p3idx, p4idx, maxDiff) {
  this.scene = scene;
  this.p1idx = p1idx;
  this.p2idx = p2idx;
  this.p3idx = p3idx;
  this.p4idx = p4idx;
  this.maxDiff = maxDiff;
  this.p12Didx = p1idx*2;
  this.p22Didx = p2idx*2;
  this.p32Didx = p3idx*2;
  this.p42Didx = p4idx*2;
  this.p13Didx = p1idx*3;
  this.p23Didx = p2idx*3;
  this.p33Didx = p3idx*3;
  this.p43Didx = p4idx*3;

  // INIT
  // compute rest length
  var p1 = scene.getPoint3D(p1idx);
  var p2 = scene.getPoint3D(p2idx);
  var p3 = scene.getPoint3D(p3idx);
  var p4 = scene.getPoint3D(p4idx);
  this.restDiff = Math.sqrt((p2[0]-p1[0])*(p2[0]-p1[0]) + (p2[1]-p1[1])*(p2[1]-p1[1]) + (p2[2]-p1[2])*(p2[2]-p1[2]));
  this.restDiff -= Math.sqrt((p4[0]-p3[0])*(p4[0]-p3[0]) + (p4[1]-p3[1])*(p4[1]-p3[1]) + (p4[2]-p3[2])*(p4[2]-p3[2]));
  this.restDiff = Math.abs(this.restDiff);
  this.restSq = this.restDiff * this.restDiff;

  // METHODS
  this.render = function() {
    var points = this.scene.points2D;
    var context = this.scene.context;

    context.beginPath();
    context.strokeStyle = '#33ee33';
    context.lineTo(points[this.p12Didx], points[this.p12Didx+1]);
    context.lineTo(points[this.p22Didx], points[this.p22Didx+1]);
    context.stroke();
    context.moveTo(points[this.p32Didx], points[this.p32Didx+1]);
    context.lineTo(points[this.p42Didx], points[this.p42Didx+1]);
    context.stroke();
  }

  this.constrain= function() {
    var points = this.scene.points3D;
    var delta1 = [points[this.p23Didx]-points[this.p13Didx], points[this.p23Didx+1]-points[this.p13Didx+1], points[this.p23Didx+2]-points[this.p13Didx+2]];
    var delta2 = [points[this.p43Didx]-points[this.p33Didx], points[this.p43Didx+1]-points[this.p33Didx+1], points[this.p43Didx+2]-points[this.p33Didx+2]];
    var delSq1 = delta1[0]*delta1[0] + delta1[1]*delta1[1] + delta1[2]*delta1[2];
    var delSq2 = delta2[0]*delta2[0] + delta2[1]*delta2[1] + delta2[2]*delta2[2];
    var delLen = Math.sqrt(delSq1) - Math.sqrt(delSq2); //TODO optimize

    if (delLen > this.maxDiff) {
      var diff = (delLen - this.maxDiff) / delLen * 0.23;
    } else if (delLen < -this.maxDiff) {
      var diff = (-delLen - this.maxDiff) / delLen * 0.23;
    } else {
      return;
    }

    for (var j=0; j<3; j++) {
      points[this.p13Didx+j] += diff * delta1[j];
      points[this.p23Didx+j] -= diff * delta1[j];
      points[this.p33Didx+j] -= diff * delta2[j];
      points[this.p43Didx+j] += diff * delta2[j];
    }
  }
}

jQuery(document).ready(function($) {
  var canvas = $('#baseCanvas');
  var scene = new Scene2D(canvas[0]);
  var state = 'none';
  var i;

  function addPoly(n) {
    //var x = Math.random() * canvas.width() - canvas.width()/2;
    //var y = Math.random() * canvas.height() - canvas.height()/2;
    var x = 0;
    var y = 0;
    var r = 50;
    n = n||4;

    var lastDir = Math.random() * 2 * Math.PI;
    var lastx = x - Math.random() * 3 * Math.cos(lastDir);
    var lasty = y - Math.random() * 3 * Math.sin(lastDir);

    var pIdx = [];
    for (var i=0; i<n; i++) {
      pIdx[i] = scene.addPoint(x+r*Math.cos(i*2*Math.PI/n), y+r*Math.sin(i*2*Math.PI/n));
    }
    scene.addItem(new polyItem(scene, pIdx));
    for (var i=0; i<n; i++) {
      scene.pointsLast[pIdx[i]*2] = lastx+r*Math.cos(i*2*Math.PI/n);
      scene.pointsLast[pIdx[i]*2+1] = lasty+r*Math.sin(i*2*Math.PI/n);
    }
  }

  //var points = {};
  //for (var p in pointData) {
  //  var point = pointData[p];
  //  points[p] = scene.addPoint(point[0], point[1]);
  //}

  //for (var item in quadItemData) {
  //  var pList = quadItemData[item]['points'];
  //  scene.addItem(new quadItem(scene, points[pList[0]], points[pList[1]], points[pList[2]], points[pList[3]]));
  //  scene.pointsLast[points[pList[0]]*2] -= quadItemData[item]['v'][0]  ;// * (0.5+Math.random());
  //  scene.pointsLast[points[pList[0]]*2+1] -= quadItemData[item]['v'][1];// * (0.5+Math.random());
  //  scene.pointsLast[points[pList[1]]*2] -= quadItemData[item]['v'][0]  ;// * (0.5+Math.random());
  //  scene.pointsLast[points[pList[1]]*2+1] -= quadItemData[item]['v'][1];// * (0.5+Math.random());
  //  scene.pointsLast[points[pList[2]]*2] -= quadItemData[item]['v'][0]  ;// * (0.5+Math.random());
  //  scene.pointsLast[points[pList[2]]*2+1] -= quadItemData[item]['v'][1];// * (0.5+Math.random());
  //  scene.pointsLast[points[pList[3]]*2] -= quadItemData[item]['v'][0]  ;// * (0.5+Math.random());
  //  scene.pointsLast[points[pList[3]]*2+1] -= quadItemData[item]['v'][1];// * (0.5+Math.random());
  //}

  //for (var item in linSumPairItemData) {
  //  var data = linSumPairItemData[item];
  //  scene.addItem(new linSumPairItem(scene, points[data[0]], points[data[1]], points[data[2]], points[data[3]], data[4], data[5]));
  //}

  //for (var item in linDiffPairItemData) {
  //  var data = linDiffPairItemData[item];
  //  scene.addItem(new linDiffPairItem(scene, points[data[0]], points[data[1]], points[data[2]], points[data[3]], data[4]));
  //}

  var frames=0;

  setInterval(function () {
    scene.tic();
    frames++;
  }, renderInterval);

  setInterval(function () {
    if (drawText) {
      framerate = frames;
      frames=0;
    }
  }, 1000);

  //setInterval(function () {
  //  if (scene.numItems < 10)
  //    addQuad();
  //}, 4000);

  $(document).keydown(function(event) {
    if (state == 'alert_next_key') {
      state = 'none';
      alert(event.which);
    }

    // KEYS
    switch(event.which) {
      case 8:                                       //backspace
        break;
      case 9:                                       //tab
        break;
      case 13:                                      //enter
        break;
      case 27:                                      //esc
        break;
      case 37:                                      //left arrow
        if (event.shiftKey) {
        } else {
          scene.panX -= 5;
        }
        break;
      case 38:                                      //up arrow
        scene.panY -= 5;
        break;
      case 39:                                      //right arrow
        if (event.shiftKey) {
        } else {
          scene.panX += 5;
        }
        break;
      case 40:                                      //down arrow
        scene.panY += 5;
        break;
      case 46:                                      //delete
        break;
      case 48:                                      //'0'
        break;
      case 56:                                      //'8' and '*'
        break;
      case 48: case 49: case 50: case 51: case 52:  //'0'-'4'
      case 53: case 54: case 55: case 56: case 57:  //'5'-'9'
        break;
      case 61:                                      //'=' / '+'
        break;
      case 66:                                      //'b' bomb
        var x = Math.random() * canvas.width() - canvas.width()/2;
        var y = Math.random() * canvas.height() - canvas.height()/2;
        var pow = 0.02;

        for (var i=0; i<scene.points3D.length; i+=3) {
          var delta = [scene.points3D[i]-x, scene.points3D[i+1]-y, scene.points3D[i+2]-z];
          var deltalength = Math.sqrt(delta[0]*delta[0] + delta[1]*delta[1] + delta[2]*delta[2]); //TODO optimize
          var c = pow/(deltalength*deltalength*deltalength);
          for (var j=0; j<3; j++) {
            scene.points3D[i+j] += delta[j]*c;
          }
        }
        break;
      case 67:                                      //'c'
        break;
      case 68:                                      //'d'
        scene.zTranslate += 0.2;
        break;
      case 69:                                      //'e'
        break;
      case 70:                                      //'f'
        addPoly(Math.floor(Math.random()*6+3));
        break;
      case 71:                                      //'g'
        gravity = !gravity;
        break;
      case 72:                                      //'h'
        hoist = !hoist;
        break;
      case 73:                                      //'i'
        scene.scale *= 1.1;
        break;
      case 76:                                      //'l'
        break;
      case 78:                                      //'n'
        drawText = !drawText;
        break;
      case 79:                                      //'o'
        scene.scale /= 1.1;
        break;
      case 80:                                      //'p'
        pause = !pause;
        break;
      case 82:                                      //'r'
        break;
      case 83:                                      //'s'
        break;
      case 85:                                      //'u'
        scene.zTranslate -= 0.2;
        break;
      case 86:                                      //'v'
        break;
      case 87:                                      //'w'
        break;
      case 89:                                      //'y'
        break;
      case 90:                                      //'z'
        break;
      case 173:                                     //'-'
        break;
      case 190:                                     //'.'
        break;
      case 191:                                     //'/'
        break;
      case 192:                                     //'`'
        if (event.shiftKey)
          state = 'alert_next_key';
        break;
      case 220:                                     //'\'
        break;
    }
  });

  // TODO fun gui stuff
  //$(document).click(function(event) {
  //});
});
