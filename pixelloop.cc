
      for (int j=0; j<npix; j++) {                          // Loop over pixels
	vector<float> skyvec;
	for (unsigned int i = 0; i < imagpx.size(); i++) {  // Loop over image stack
	  float pval = imagpx[i][j];                        // Data value
	  int mval = maskpx[i][j];                          // Mask value 
	  
	  // Insert into a sorted list unless the pixel is masked
	  //
	  if (mval == 0) {                                  // Not masked
	    bool inserted = false;
	    vector<float>::iterator ii = skyvec.begin();
	    while (!inserted && skyvec.size() > 0 && ii != skyvec.end()) {
	      if (pval < *ii) {
		skyvec.insert(ii,pval);
		inserted = true;
	      }
	      ii++;
	    }
	    if (!inserted) skyvec.push_back(pval);          // Insert at the end
	  }
	}

	int nim = skyvec.size() - 2*nrej;   // Number of images used for this pixel
	if (nim < 0) nim = 0;
	nimuse[j] = nim;

	// Are there enough values in the list to continue?
	//
	if (nim > 1) {
	  int n1=nrej, n2=skyvec.size()-nrej;

	  // Calculate mean
	  //
	  float val = 0.0;
	  for (int n=n1; n<n2; n++) {
	    float foo = skyvec[n];
	    val += foo;
	  }
	  val /= n2-n1;
	  skyval[j] = val;

	  // Calculate RMS
	  //
	  float rms = 0.0;
	  for (int n=n1; n<n2; n++) {
	    float foo = skyvec[n];
	    foo -= val;
	    rms += foo * foo;
	  }
	  rms /= n2-n1;
	  skyrms[j] = rms;
	} else {
	  skyval[j] = 0.0;
	  skyrms[j] = 0.0;
	}

      } // End of loop over pixels
