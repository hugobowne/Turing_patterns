#include "ofApp.h"

//--------------------------------------------------------------
void ofApp::setup(){

	ofSetVerticalSync(true);
	ofBackground(20);
    
    // GL_REPEAT for texture wrap only works with NON-ARB textures //
    ofDisableArbTex();
    texture.loadImage("of.png");
    texture.getTextureReference().setTextureWrap( GL_REPEAT, GL_REPEAT );
    vidGrabber.initGrabber(640, 480, true);
    
    bFill       = true;
    bWireframe  = true;
    bDrawNormals= false;
    bDrawAxes   = false;
    bDrawLights = false;
    bInfoText   = true;
    bMousePressed   = false;
    bSplitFaces = false;
    
    float width     = ofGetWidth() * .2;
    float height    = ofGetHeight() * .2;
    

    icoSphere.setRadius( width );
    
    mode = 0;
    
    ofSetSmoothLighting(true);
    pointLight.setDiffuseColor( ofFloatColor(.85, .85, .55) );
    pointLight.setSpecularColor( ofFloatColor(1.f, 1.f, 1.f));
    
    pointLight2.setDiffuseColor( ofFloatColor( 238.f/255.f, 57.f/255.f, 135.f/255.f ));
    pointLight2.setSpecularColor(ofFloatColor(.8f, .8f, .9f));
    
    pointLight3.setDiffuseColor( ofFloatColor(19.f/255.f,94.f/255.f,77.f/255.f) );
    pointLight3.setSpecularColor( ofFloatColor(18.f/255.f,150.f/255.f,135.f/255.f) );
    
    // shininess is a value between 0 - 128, 128 being the most shiny //
	material.setShininess( 120 );
    // the light highlight of the material //
	material.setSpecularColor(ofColor(255, 255, 255, 255));
    
    ofSetSphereResolution(24);
    
}

//--------------------------------------------------------------
void ofApp::update() {
    pointLight.setPosition((ofGetWidth()*.5)+ cos(ofGetElapsedTimef()*.5)*(ofGetWidth()*.3), ofGetHeight()/2, 500);
    pointLight2.setPosition((ofGetWidth()*.5)+ cos(ofGetElapsedTimef()*.15)*(ofGetWidth()*.3),
                            ofGetHeight()*.5 + sin(ofGetElapsedTimef()*.7)*(ofGetHeight()), -300);
    
    pointLight3.setPosition(
                            cos(ofGetElapsedTimef()*1.5) * ofGetWidth()*.5,
                            sin(ofGetElapsedTimef()*1.5f) * ofGetWidth()*.5,
                            cos(ofGetElapsedTimef()*.2) * ofGetWidth()
    );
    
	//ofSetWindowTitle("Framerate: "+ofToString(ofGetFrameRate(), 0));
    if(mode == 2 || ofGetElapsedTimef() < 10) {
        vidGrabber.update();
    }
}

//--------------------------------------------------------------
void ofApp::draw() {
    
    float spinX = sin(ofGetElapsedTimef()*.35f);
    float spinY = cos(ofGetElapsedTimef()*.075f);
    
    if(bMousePressed) {
        spinX = spinY = 0.0f;
    }
    
    ofEnableDepthTest();
    
    ofEnableLighting();
    pointLight.enable();
    pointLight2.enable();
    pointLight3.enable();
    
	material.begin();
    
    
    ofSetColor(180);
    ofNoFill();
    ofDrawSphere(ofGetWidth()/2, ofGetHeight()/2, ofGetWidth());
    
    if(mode == 1 || mode == 3) texture.getTextureReference().bind();
    if(mode == 2) vidGrabber.getTextureReference().bind();
    
    
    
    
    // ICO Sphere //
    icoSphere.setPosition(ofGetWidth()*.5, ofGetHeight()*.5, 0);
    icoSphere.rotate(spinX, 1.0, 0.0, 0.0);
    icoSphere.rotate(spinY, 0, 1.0, 0.0);
    
    if(mode == 3) {
        triangles = icoSphere.getMesh().getUniqueFaces();
    }
    
    if(bFill) {
        ofFill();
        ofSetColor(255);
        
        if(mode == 3) {
            float angle = (ofGetElapsedTimef() * 1.4);
            ofVec3f faceNormal;
            for(int i = 0; i < triangles.size(); i++ ) {
                float frc = ofSignedNoise(angle* (float)i * .1, angle*.05) * 4;
                faceNormal = triangles[i].getFaceNormal();
                for(int j = 0; j < 3; j++ ) {
                    triangles[i].setVertex(j, triangles[i].getVertex(j) + faceNormal * frc );
                }
            }
            icoSphere.getMesh().setFromTriangles( triangles );
        }
        
        icoSphere.draw();
    }
    if(bWireframe) {
        ofNoFill();
        ofSetColor(0, 0, 0);
        if(!bFill) ofSetColor(255);
        icoSphere.setScale(1.01f);
        icoSphere.drawWireframe();
        icoSphere.setScale(1.f);
    }

    
    if(mode == 1 || mode == 3) texture.getTextureReference().unbind();
    if(mode == 2) vidGrabber.getTextureReference().unbind();
    
    material.end();
    ofDisableLighting();
    
    if(bDrawLights) {
        ofFill();
        ofSetColor(pointLight.getDiffuseColor());
        pointLight.draw();
        ofSetColor(pointLight2.getDiffuseColor());
        pointLight2.draw();
        ofSetColor(pointLight3.getDiffuseColor());
        pointLight3.draw();
    }
    
    if(bDrawNormals) {
        ofSetColor(225, 0, 255);
        icoSphere.drawNormals(20, bSplitFaces);
    }
    if(bDrawAxes) {
        icoSphere.drawAxes(icoSphere.getRadius()+30);
    }
    
    ofDisableDepthTest();
    
    ofFill();
    
    ofSetColor(0);
    ofRect(icoSphere.getPosition().x-154, icoSphere.getPosition().y + 120, 168, 24);
    ofSetColor(255);
    ofDrawBitmapString("ofIcoSpherePrimitive", icoSphere.getPosition().x-150, icoSphere.getPosition().y+136 );
    
    if(bInfoText) {
        stringstream ss;
        ss << "Framerate: " << ofToString(ofGetFrameRate(),0) << "\n";
        ss << "(f): Toggle Fullscreen"<<endl<<"(s): Draw Solid Shapes"<<endl<<"(w): Draw Wireframes"<<endl;
        ss <<"(1/2/3/4): Set Resolutions" <<endl<<"(n): Draw Normals"<<"\n(LEFT/RIGHT): Set Mode "<<ofToString(mode,0)<<endl;
        ss <<"(z): Split Faces " <<bSplitFaces<<endl;
        ss <<"(a): Draw Axes"<<endl<<"(l): Render lights"<<endl<<"(t): Info Text"<<endl;
        
        ofDrawBitmapString(ss.str().c_str(), 20, 20);
    }
    
    
}

//--------------------------------------------------------------
void ofApp::keyPressed(int key) {
    
	switch(key) {
		case 'f':
			ofToggleFullscreen();
			break;
        case 's':
            bFill = !bFill;
            break;
        case 'w':
            bWireframe = !bWireframe;
            break;
        case '1':
            bSplitFaces=false;
            // icosahedron //
            icoSphere.setResolution(0); // number of subdivides //
            break;
        case '2':
            bSplitFaces=false;
            icoSphere.setResolution(1);
            break;
        case '3':
            bSplitFaces=false;
            icoSphere.setResolution(2);
            break;
        case '4':
            bSplitFaces=false;
            icoSphere.setResolution(4);
            break;
        case 'n':
            bDrawNormals = !bDrawNormals;
            break;
        case OF_KEY_RIGHT:
            mode++;
            if(mode > 3) mode = 0;
            if(mode==3){
                // to get unique triangles, you have to use triangles mode //
                //sphere.setMode( OF_PRIMITIVE_TRIANGLES );
            }
            break;
        case OF_KEY_LEFT:
            mode--;
            if(mode < 0) mode = 3;
            if(mode==3){
                // to get unique triangles, you have to use triangles mode //
                //sphere.setMode( OF_PRIMITIVE_TRIANGLES );
            }
            break;
        case 'a':
            bDrawAxes = !bDrawAxes;
            break;
        case 'l':
            bDrawLights = !bDrawLights;
            break;
        case 't':
            bInfoText=!bInfoText;
            break;
        case 'z':
            bSplitFaces = !bSplitFaces;
            
            if(mode == 3) bSplitFaces = false;
            
            if(bSplitFaces) {
                icoSphere.getMesh().setFromTriangles(triangles, true);
                
            } else {
                // vertex normals are calculated with creation, set resolution //
                
                icoSphere.setResolution( icoSphere.getResolution() );
            }
            break;
	}
    
    if(mode == 1) {
        // resize the plane to the size of the texture //
        // setTexCoordsFromTexture sets normalized or non-normalized tex coords based on an ofTexture passed in.
        icoSphere.mapTexCoordsFromTexture( texture.getTextureReference() );
    }
    
    if(mode == 2) {
        icoSphere.mapTexCoordsFromTexture( vidGrabber.getTextureReference() );
    }
    
    // 
    if( mode == 3 ) {
        
        bSplitFaces = false;
        
        // if the faces were split, we can get some weird results, since we
        // might not know what the new strides were,
        // so reset the primitives by calling their setMode function
        // which recreates the mesh with the proper indicies //
        icoSphere.setMode(OF_PRIMITIVE_TRIANGLE_STRIP);
        // box only supports triangles //
        
        
        // rebuild the box, 
        icoSphere.mapTexCoords(0, 0, 5, 5);
        
    }

}

//--------------------------------------------------------------
void ofApp::keyReleased(int key){

}

//--------------------------------------------------------------
void ofApp::mouseMoved(int x, int y ){
	
}

//--------------------------------------------------------------
void ofApp::mouseDragged(int x, int y, int button){


}

//--------------------------------------------------------------
void ofApp::mousePressed(int x, int y, int button){
    bMousePressed = true;
}

//--------------------------------------------------------------
void ofApp::mouseReleased(int x, int y, int button){
    bMousePressed = false;
}

//--------------------------------------------------------------
void ofApp::windowResized(int w, int h){

}

//--------------------------------------------------------------
void ofApp::gotMessage(ofMessage msg){

}

//--------------------------------------------------------------
void ofApp::dragEvent(ofDragInfo dragInfo){ 

}
