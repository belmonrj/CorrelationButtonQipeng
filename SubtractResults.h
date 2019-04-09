class subResult {

  public:
    subResult() {
        m_coeff_raw_value.clear();
        m_coeff_raw_error.clear();
        m_coeff_sub_value.clear();
        m_coeff_sub_error.clear();
        m_coeff_subImp_value.clear();
        m_coeff_subImp_error.clear();
    };
    virtual ~subResult() {};

    void setCoeffRaw(vector<float>& vec_value, vector<float>& vec_error) { 
        m_coeff_raw_value = vec_value; 
        m_coeff_raw_error = vec_error; 
    }

    void setCoeffSub(vector<float>& vec_value, vector<float>& vec_error) { 
        m_coeff_sub_value = vec_value; 
        m_coeff_sub_error = vec_error; 
    }

    void setCoeffSubMinos(vector<float>& vec_error_lower, vector<float>& vec_error_upper) { 
        m_coeff_sub_lowerError = vec_error_lower; 
        m_coeff_sub_upperError = vec_error_upper; 
    }

    void setCoeffSubImp(vector<float>& vec_value, vector<float>& vec_error) { 
        m_coeff_subImp_value = vec_value; 
        m_coeff_subImp_error = vec_error; 
    }

    // methods
    // most common methods for just v22
    float getV22RawValue() { if (m_coeff_raw_value.size() > 1) return m_coeff_raw_value[1]; else return 0.;}
    float getV22RawError() { if (m_coeff_raw_error.size() > 1) return m_coeff_raw_error[1]; else return 0.;}
    float getV22SubValue() { if (m_coeff_sub_value.size() > 1) return m_coeff_sub_value[1]; else return 0.;}
    float getV22SubError() { if (m_coeff_sub_error.size() > 1) return m_coeff_sub_error[1]; else return 0.;}

    float getV22SubImpValue() { if (m_coeff_subImp_value.size() > 1) return m_coeff_subImp_value[1]; else return 0.;}
    float getV22SubImpError() { if (m_coeff_subImp_error.size() > 1) return m_coeff_subImp_error[1]; else return 0.;}

    float getChi2() {return m_fitChi2;}
    void  setChi2(float _chi2) {m_fitChi2 = _chi2;}

    // advanced usage
    void  setNHar(int _n) {m_nhar = _n;}
    int   getNHar()     {return m_nhar;}
    float getRhoValue() {return m_rho_value;}
    float getRhoError() {return m_rho_error;}
    void  setRhoValue(float _val) {m_rho_value = _val;}
    void  setRhoError(float _err) {m_rho_error = _err;}

    float getPedstalValue() {return m_pedstalG;}
    float getPedstalError() {return m_pedstalG_error;}
    void  setPedstalValue(float _val) {m_pedstalG = _val;}
    void  setPedstalError(float _err) {m_pedstalG_error = _err;}

    float getCoeffRawValue(int order) {
        if (order < m_nhar+1 && m_coeff_raw_value.size() > order-1) return m_coeff_raw_value[order-1];
        else return 0;
    }

    float getCoeffRawError(int order) {
        if (order < m_nhar+1 && m_coeff_raw_value.size() > order-1) return m_coeff_raw_error[order-1];
        else return 0;
    }

    float getCoeffSubValue(int order) {
        if (order < m_nhar+1 && m_coeff_sub_value.size() > order-1) return m_coeff_sub_value[order-1];
        else return 0;
    }

    float getCoeffSubError(int order) {
        if (order < m_nhar+1 && m_coeff_sub_value.size() > order-1) return m_coeff_sub_error[order-1];
        else return 0;
    }

    float getCoeffSubLowerError(int order) {
        if (order < m_nhar+1 && m_coeff_sub_lowerError.size() > order-1) return m_coeff_sub_lowerError[order-1];
        else return 0;
    }

    float getCoeffSubUpperError(int order) {
        if (order < m_nhar+1 && m_coeff_sub_upperError.size() > order-1) return m_coeff_sub_upperError[order-1];
        else return 0;
    }

    float getCoeffSubImpValue(int order) {
        if (order < m_nhar+1 && m_coeff_subImp_value.size() > order-1) return m_coeff_subImp_value[order-1];
        else return 0;
    }

    float getCoeffSubImpError(int order) {
        if (order < m_nhar+1 && m_coeff_subImp_value.size() > order-1) return m_coeff_subImp_error[order-1];
        else return 0;
    }

    float getRhoCorrelation(int order) {
        if (order < m_nhar+1 && m_correlation_rho_flowCoeff.size() > order-1) return m_correlation_rho_flowCoeff[order-1];
        else return 0;
    }

    void setRhoCorrelation(vector<float>& vec) { 
        m_correlation_rho_flowCoeff = vec; 
    }
    
  private:
    // members of flow coefficients
    // three types to cover usage of different methods
    // set to zero if certain type is not availble

    // unsubtracted flow coefficients
    vector<float> m_coeff_raw_value;
    vector<float> m_coeff_raw_error;

    // subtracted flow coefficients
    vector<float> m_coeff_sub_value;
    vector<float> m_coeff_sub_error;
    // introduce members for handling Minos errors
    // not sure if it's needed for raw or subImp coefficients
    vector<float> m_coeff_sub_lowerError;
    vector<float> m_coeff_sub_upperError;

    // subtracted and corrected (for LM flow) flow coefficeints
    vector<float> m_coeff_subImp_value;
    vector<float> m_coeff_subImp_error;

    int m_nhar;

    // additional common fitting parameters
    float m_pedstalG; //G_HM in case of ATLAS fit
    float m_pedstalG_error;

    float m_fitChi2;
    float m_rho_value;
    float m_rho_error;
    vector<float> m_correlation_rho_flowCoeff;
};
