import numpy as np
def ukca_ht():
	height =[]
	z_top_of_model=85000.00 # 85km zT
	h=100.0
	firstconstantrrholevel=51,
	# eta_theta = n; eta_rho = nI
	eta_theta=np.array([0.000000e+00, 0.2352941e-03,0.6274510e-03,0.1176471e-02,0.1882353e-02,
	0.2745098e-02,0.3764706e-02,0.4941176e-02,0.6274510e-02,0.7764705e-02,0.9411764e-02,0.1121569e-01,0.1317647e-01,0.1529412e-01,0.1756863e-01,
	0.2000000e-01,0.2258823e-01,0.2533333e-01,0.2823529e-01,0.3129411e-01,0.3450980e-01,0.3788235e-01,0.4141176e-01,0.4509804e-01,0.4894118e-01,
	0.5294117e-01,0.5709804e-01,0.6141176e-01,0.6588235e-01,0.7050980e-01,0.7529411e-01,0.8023529e-01,0.8533333e-01,0.9058823e-01,0.9600001e-01,
	0.1015687e+00,0.1072942e+00,0.1131767e+00,0.1192161e+00,0.1254127e+00,0.1317666e+00,0.1382781e+00,0.1449476e+00,0.1517757e+00,0.1587633e+00,
	0.1659115e+00,0.1732221e+00,0.1806969e+00,0.1883390e+00,0.1961518e+00,0.2041400e+00,0.2123093e+00,0.2206671e+00,0.2292222e+00,0.2379856e+00,
	0.2469709e+00,0.2561942e+00,0.2656752e+00,0.2754372e+00,0.2855080e+00,0.2959203e+00,0.3067128e+00,0.3179307e+00,0.3296266e+00,0.3418615e+00,
	0.3547061e+00,0.3682416e+00,0.3825613e+00,0.3977717e+00,0.4139944e+00,0.4313675e+00,0.4500474e+00,0.4702109e+00,0.4920571e+00,0.5158098e+00,
	0.5417201e+00,0.5700686e+00,0.6011688e+00,0.6353697e+00,0.6730590e+00,0.7146671e+00,0.7606701e+00,0.8115944e+00,0.8680208e+00,0.9305884e+00, 0.10000e+01])

	eta_rho=np.array([0.1176471e-03,0.4313726e-03,0.9019608e-03,0.1529412e-02,0.2313725e-02,0.3254902e-02,0.4352941e-02,0.5607843e-02,0.7019607e-02,
	0.8588235e-02,0.1031373e-01,0.1219608e-01,0.1423529e-01,0.1643137e-01,0.1878431e-01,0.2129412e-01,0.2396078e-01,0.2678431e-01,0.2976470e-01,
	0.3290196e-01,0.3619608e-01,0.3964706e-01,0.4325490e-01,0.4701960e-01,0.5094118e-01,0.5501961e-01,0.5925490e-01,0.6364705e-01,0.6819607e-01,
	0.7290196e-01,0.7776470e-01,0.8278431e-01,0.8796078e-01,0.9329412e-01,0.9878433e-01,0.1044314e+00,0.1102354e+00,0.1161964e+00,0.1223144e+00,
	0.1285897e+00,0.1350224e+00,0.1416128e+00,0.1483616e+00,0.1552695e+00,0.1623374e+00,0.1695668e+00,0.1769595e+00,0.1845180e+00,0.1922454e+00,
	0.2001459e+00,0.2082247e+00,0.2164882e+00,0.2249446e+00,0.2336039e+00,0.2424783e+00,0.2515826e+00,0.2609347e+00,0.2705562e+00,0.2804726e+00,
	0.2907141e+00,0.3013166e+00,0.3123218e+00,0.3237787e+00,0.3357441e+00,0.3482838e+00,0.3614739e+00,0.3754014e+00,0.3901665e+00,0.4058831e+00,
	0.4226810e+00,0.4407075e+00,0.4601292e+00,0.4811340e+00,0.5039334e+00,0.5287649e+00,30.5558944e+00,0.5856187e+00,0.6182693e+00,0.6542144e+00,
	0.6938630e+00,0.7376686e+00,0.7861323e+00,0.8398075e+00,0.8993046e+00,0.9652942e+00])

	for i in range (1, 86):
		if eta_theta[i] <= eta_rho[firstconstantrrholevel]:
			height.append((eta_theta[i]*z_top_of_model)+ (h*(1-((eta_theta[i]/eta_rho[firstconstantrrholevel])**2))))
		else:
			height.append(eta_theta[i]*z_top_of_model)
	# print len(height)
	height_a = height[:52]
	print len(height_a)
	return height_a
	
 
def height_id():
	id = []
	altitude = ukca_ht()
	
	for idx in range(0, 51):
		if 0<altitude[idx]<=altitude[idx+1]:
			id.append(idx+1)
		# print id
	return id
	
print(ukca_ht())
print(height_id())