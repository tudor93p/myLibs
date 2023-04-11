using myLibs: Utils
using BenchmarkTools

import PyPlot  



#Utils.Distribute_Work(1:3,(args...;kwargs...)->sleep(2*rand()))

@testset "Utils.reduce_index=mod1" begin 



for n=1:100,i=-300:300 

#   @test Utils.reduce_index(i,n)==mod1(i,n)

end
end 


 
@testset "path connect" begin 

	points0 = sort(rand(3))

	points0 = [0,1,7]


	dim = 2 
	
	points = reshape(points0,1,length(points0))  


	path, xticks, x = Utils.PathConnect(points, 13; dim=dim, bounds=[1.2,5])

	PyPlot.scatter(x, vcat(path...))





end 





@testset "periodic distance types" begin 

	for a in (rand(1:5), rand(),rand(3),rand(1:5,3) ),T in (rand(),6)

		for b in (rand(1:5), rand(),rand(3),rand(1:5,3) )

			@test applicable(Utils.dist_periodic, a, b, T)
			
	
#			@show a b T 

			for (d1,d2,d3,d4) in zip(
														Utils.dist_periodic(a, b, T),
														Utils.dist_periodic(Float64.(a), b, T),
														Utils.dist_periodic(a, Float64.(b), T),
														Utils.dist_periodic(a, b, Float64(T)),
														Float64.(Utils.dist_periodic(a, b, T)),
														)

	@test d1≈d2 ≈d3≈d4

			end 



		end 

	end 

end 



@testset "periodic distance" begin 

	@test Utils.dist_periodic([1,2,3],1,3)==[0,1,1]

	@test Utils.dist_periodic(1,2,1)==0
	@test Utils.dist_periodic(1,1.1,1)≈0.1
	@test Utils.dist_periodic(1.1pi,3pi,2pi)≈0.1pi

	A = rand(2,3) .- 0.5 
	B = rand(2,3) .- 0.5 

	@test Utils.dist_periodic(A,B,0.1,)≈Utils.dist_periodic(B,A,0.1,)

	@test Utils.dist_periodic(A[1],B,0.1,)≈Utils.dist_periodic(B,A[1],0.1,)
	@test Utils.dist_periodic(A[1],B,0.1,)[1]≈Utils.dist_periodic(B,A,0.1,)[1]


end 


error()

@testset "closest data points" begin 

	X = sort(rand(100))

	for i in 5:95 

		i1,i2 = Utils.closest_data_points(X, X[i]/2+X[i+1]/2, 3)

		@test i1 == i-2:i
		@test i2 == i+1:i+3

	end 


	@show Utils.closest_data_points(X, 1.5, 3)
	@show Utils.closest_data_points(X, X[50], 3)

	@show Utils.closest_data_points(sort(rand(6)), 0.5, 10)
	@show Utils.closest_data_points(rand(6), 0.5, 10)

end 




x = [0.8405599073149947, 0.8615739049978697, 0.8825879026807445, 0.9036019003636193, 0.9246158980464944, 0.9456298957293693, 0.9666438934122441, 0.9876578910951189, 1.0086718887779937, 1.0296858864608687, 1.0506998841437436, 1.0717138818266183, 1.0927278795094932, 1.1137418771923682, 1.1347558748752429, 1.1557698725581178, 1.1767838702409927, 1.1977978679238677, 1.470979837801241, 1.4919938354841158, 1.5130078331669907, 1.5340218308498657, 1.5550358285327404, 1.5760498262156153, 1.5970638238984902, 1.618077821581365, 1.6390918192642399, 1.6601058169471148, 1.6811198146299895, 1.7021338123128644, 1.7231478099957394, 1.744161807678614, 1.765175805361489, 1.786189803044364, 1.8072038007272386, 1.8282177984101136, 1.8492317960929887, 1.8702457937758636, 1.8912597914587386, 1.9122737891416133, 1.9332877868244882, 1.9543017845073631, 1.9753157821902378, 1.9963297798731128, 2.0173437775559875, 2.0383577752388624, 2.0593717729217373, 2.0803857706046123, 2.101399768287487, 2.1224137659703617, 2.1434277636532366, 2.1644417613361115, 2.1854557590189865, 2.2064697567018614, 2.2274837543847363, 2.248497752067611, 2.2695117497504858, 2.2905257474333607, 2.3115397451162356, 2.3325537427991105, 2.3535677404819855, 2.3745817381648604, 2.3955957358477353, 2.4166097335306103, 2.437623731213485, 2.45863772889636, 2.4796517265792346, 2.5006657242621095, 2.5216797219449845, 2.5426937196278594, 2.5637077173107343, 2.5847217149936093, 2.6057357126764837, 2.6267497103593587, 2.6477637080422336, 2.6687777057251085, 2.6897917034079835, 2.7108057010908584, 2.731819698773733, 2.752833696456608, 2.7738476941394827, 2.7948616918223577, 2.8158756895052326, 2.8368896871881075, 2.857903684870982, 2.878917682553857, 2.8999316802367323, 2.9209456779196072, 2.941959675602482, 2.9629736732853567, 2.9839876709682316, 3.0050016686511065, 3.0260156663339814, 3.0470296640168564, 3.0680436616997313, 3.089057659382606, 3.1100716570654807, 3.1310856547483557, 3.1520996524312306, 3.1731136501141055, 3.1941276477969804, 3.215141645479855, 3.23615564316273, 3.257169640845605, 3.2781836385284797, 3.2991976362113546, 3.3202116338942296, 3.341225631577104, 3.362239629259979, 3.383253626942854, 3.404267624625729, 3.425281622308604, 3.4462956199914787, 3.467309617674353, 3.488323615357228, 3.509337613040103, 3.530351610722978, 3.551365608405853, 3.572379606088728, 3.5933936037716023, 3.6144076014544773, 3.635421599137352, 3.656435596820227, 3.6774495945031025, 3.6984635921859774, 3.7194775898688524, 3.7404915875517273, 3.761505585234602, 3.782519582917477, 3.8035335806003516, 3.8245475782832266, 3.8455615759661015, 3.8665755736489764, 3.8875895713318513, 3.9086035690147263, 3.9296175666976008, 3.9506315643804757, 3.9716455620633506, 3.9926595597462256, 4.0136735574291, 4.034687555111975, 4.05570155279485, 4.076715550477725, 4.0977295481606, 4.118743545843475, 4.13975754352635, 4.1607715412092245, 4.1817855388920995, 4.202799536574974, 4.223813534257849, 4.244827531940723, 4.265841529623598, 4.286855527306473, 4.307869524989348, 4.328883522672223, 4.349897520355098, 4.370911518037973, 4.391925515720848, 4.412939513403723, 4.433953511086598, 4.454967508769473, 4.475981506452348, 4.496995504135222, 4.518009501818097, 4.5390234995009715, 4.560037497183846, 4.581051494866721, 4.602065492549596, 4.623079490232471, 4.644093487915346, 4.665107485598221, 4.686121483281096, 4.707135480963971, 4.728149478646846];

y = [-0.28553738104810655, -0.25517429525520985, -0.22442256003846506, -0.1933638565729323, -0.16212995843793632, -0.1309660233032135, -0.10042240012629389, -0.07208388865104709, -0.051853817204905894, 0.032240262261239216, 0.0547727536879421, 0.08438506432738357, 0.11614962833689013, 0.14862698538550934, 0.18118585343886906, 0.21338818011532545, 0.24479604649179798, 0.27486872122833805, -0.27519724426979764, -0.24587776940822614, -0.21884852758187012, -0.1948843730146601, -0.28774542102218437, -0.2970141292452641, -0.29990245192429715, -0.29967600338983874, -0.29929916840439685, -0.2988093312042299, -0.298177340886885, -0.2974105185642764, -0.2965239965251948, -0.29550591154622774, -0.2943506289142895, -0.29307731326830416, -0.29165760679414143, -0.2901332031942324, -0.2884751422878602, -0.28668937044123055, -0.2847212912490798, -0.2826796354946793, -0.28050656489384435, -0.2782024605328075, -0.2757905486311976, -0.27324958153003964, -0.27057861409931494, -0.26780342434729904, -0.2649085161087776, -0.2618838410708863, -0.2587480947819058, -0.2555130592910209, -0.2521540262551555, -0.24867093997952086, -0.24509172899172138, -0.2414173707514387, -0.23762260189054155, -0.23370720041363888, -0.22970018994916488, -0.22561491871530212, -0.2214225900724216, -0.21710715789580398, -0.21268761712262665, -0.20820679878660117, -0.2036633854440101, -0.19902109281069352, -0.19424636453538424, -0.18935015032221855, -0.1844559807991386, -0.18204153229205255, -0.18159018243635483, -0.18117134699878967, -0.18077360438466838, -0.18043520772227994, -0.18018054200737676, -0.1800108727526446, -0.17991355329248293, -0.1798764996707769, -0.1309770306051463, -0.1252758349765525, -0.11951849476492761, -0.11370915609082019, -0.10785254834357072, -0.10194542234590206, -0.09599410788153098, -0.09000234418910254, -0.08397084695652744, -0.07790186596766523, -0.07179792005630152, -0.06566171826860236, -0.059496090043243, -0.053303918157349965, -0.047088089500887435, -0.040851467850178416, -0.03459688500734444, -0.028327143732073033, -0.02204502625450942, -0.015753303970393092, -0.009454745941372782, -0.0031521260655557513, 0.003151769922426086, 0.009454148027606486, 0.01575220523069584, 0.02204312477809262, 0.028324083398219388, 0.03459226144504486, 0.04084485908342669, 0.04707912058790372, 0.053292365862270825, 0.05948202365672504, 0.06564565383210026, 0.0717809377098072, 0.0778856101618705, 0.08395731197791802, 0.08999336384324326, 0.095990478524805, 0.10194331973143414, 0.10784907066406565, 0.11370737264111778, 0.11950895250633477, 0.12525417609398332, 0.13094505509333002, 0.13658434368759848, 0.1421728348440858, 0.1477052211696569, 0.15316564994868084, 0.15854636672946368, 0.16382487981641683, 0.1690048756050031, 0.17410330880595473, 0.17909717956687907, 0.1820537869355426, 0.18246011395813122, 0.1827772777598664, 0.18309919212209547, 0.18341382476472168, 0.18368047653183384, 0.18386663044228557, 0.1839667255138333, 0.22138848306858636, 0.22559017570734366, 0.22969205449454966, 0.2336964678935417, 0.2376023797822959, 0.24139764136373273, 0.24508080312868002, 0.2486617496340141, 0.2521391450273379, 0.25549811596871813, 0.25873765883733985, 0.26187182008088183, 0.26489747261559954, 0.2677892702459472, 0.2705651221736405, 0.27323811968110184, 0.275776742115873, 0.27817906666651354, 0.2804903519011911, 0.2826571442192366, 0.284605429001025, 0.28664938246960414, 0.28844391290855803, 0.29011066897845533, 0.2916518799826605, 0.29306068957991177, 0.2943444742201416, 0.2954948645869532, 0.2965157075581023, 0.2973993486197317, 0.2981736494393564, 0.29878959979806613, -0.2998746347291381, -0.2949401591864813, -0.2931977286734355, 0.15994134360062295, 0.17489399333248867];
z = [-0.9922872214419209, -0.9899173058797451, -0.9863710496432058, -0.9807408969261279, -0.9710648437847998, -0.9524897342724126, -0.9105794656039895, -0.7906789861787584, -0.36412543420043547, -0.4473924787077033, -0.8105428353127486, -0.9146062435471397, -0.9527402039918733, -0.9702684089427724, -0.979659624699121, -0.9852511452682069, -0.9888411150562566, -0.9912747897244909, 0.9764628775449535, 0.9681735644155441, 0.9555858491156365, 0.9354108566131017, -0.9455493391217877, 0.9324563498142809, 0.9234067892619998, 0.921121091225792, -0.9243225608629483, -0.9134440585119096, -0.9500470276763349, -0.496325458469792, 0.9801574101062906, -0.9869084238864444, 0.9876667159749162, -0.9876361931663519, 0.9452591585134715, -0.9619118430405087, -0.8842487512107406, -0.060114912688130115, 0.986948402162839, 0.6820401687565134, 0.7161457125170307, 0.9872762771272063, 0.7672422720347641, 0.8084594327500251, 0.9716255210341608, 0.7971099867965628, 0.9793950736629703, 0.8137199674415455, 0.8836659831784154, 0.8947013829662204, 0.923578789848206, 0.8246041093047777, 0.7875554583014065, 0.9337690997262633, 0.9511609110527076, 0.7529400988252614, 0.7581302880706785, 0.8936014780620531, 0.9873910973636959, 0.8806419722448812, 0.6404308254112475, 0.6991553874301464, 0.9403854308105938, 0.9869238684456778, 0.9251762142393848, 0.44608525462739335, 0.20754365758262203, 0.9250341949153712, 0.933141057514106, 0.9358562078219771, 0.9387142469833769, 0.941678519434733, 0.9445437223571235, 0.947207077951124, 0.9496783551787091, 0.9520170211387976, -0.9875561154667388, -0.9796490638569664, -0.9386643363489092, 0.013178726694566028, 0.9867447647802557, 0.8191909614619305, -0.849165234527923, -0.9654814699219717, -0.9816108777929692, -0.9861131938283627, -0.9874852440395613, -0.9876756013965715, -0.9873278140816756, -0.9866950944484473, -0.9858772229021114, -0.9848933645384215, -0.9836998641760939, -0.9821844312268222, -0.9801763672715279, -0.9777125805488177, -0.9768685754839214, -0.9845269307746437, -0.9845268124009782, -0.9768682384825401, -0.9777120819312823, -0.9801758100807332, -0.98218397348316, -0.983699724531674, -0.9848938166113399, -0.9858785676103953, -0.9866975992048719, -0.9873316126468187, -0.9876805623656781, -0.9874908416631648, -0.9861184486467716, -0.9816144615989996, -0.9654820605458466, -0.8491621234399372, 0.8191971822195505, 0.9867518371972572, 0.013183205906087915, -0.9386654805123273, -0.9796562660377512, -0.9875661971662878, -0.9832158205185636, -0.953157747234566, -0.7196342085908537, 0.9528191752120003, -0.46253364305463107, -0.94413961198102, -0.9856865525539059, -0.976399375643215, -0.8907361115734872, -0.9023687542692402, -0.9189455251803574, -0.917707665156646, -0.9140027791633635, -0.9093570336721858, -0.9045587126783967, -0.899959800754287, -0.8954201359663616, 0.9875912540064845, 0.8936393189098162, 0.7579389130443069, 0.7528157432848147, 0.9512680377393854, 0.9338669266150489, 0.7874317810157164, 0.8244736682211816, 0.9236716237136126, 0.8947853068305565, 0.883517258355811, 0.8136358972462998, 0.9795456795222028, 0.7970836227880662, 0.9714090393186837, 0.8085824190175359, 0.7673450167960026, 0.9868553366492934, 0.7163420460925801, 0.6824496027139999, 0.9811466334554296, -0.050138763717907944, -0.8844794075408421, -0.9620337397822479, 0.9452964646010986, -0.9876630729172662, 0.9876724717138701, -0.9869080890723293, 0.9801405074903653, -0.49549071791301796, -0.9501391137738702, -0.9134710877120601, 0.41251854419535516, 0.8516008201503393, 0.9189899818126309, -0.8366022243113383, -0.9008385399132073];


x2 = [1.5550358285327404, 1.5760498262156153, 1.5970638238984902, 1.618077821581365, 1.6390918192642399, 1.6601058169471148, 1.6811198146299895, 1.7021338123128644, 1.7231478099957394, 1.744161807678614, 1.765175805361489, 1.786189803044364, 1.8072038007272386, 1.8282177984101136, 1.8492317960929887, 1.8702457937758636, 1.8912597914587386, 1.9122737891416133, 1.9332877868244882, 1.9543017845073631, 1.9753157821902378, 1.9963297798731128, 2.0173437775559875, 2.0383577752388624, 2.0593717729217373, 2.0803857706046123, 2.101399768287487, 2.1224137659703617, 2.1434277636532366, 2.1644417613361115, 2.1854557590189865, 2.2064697567018614, 2.2274837543847363, 2.248497752067611, 2.2695117497504858, 2.2905257474333607, 2.3115397451162356, 2.3325537427991105, 2.3535677404819855, 2.3745817381648604, 2.3955957358477353, 2.4166097335306103, 2.437623731213485, 2.45863772889636, 2.4796517265792346, 2.5006657242621095, 2.5216797219449845, 2.5426937196278594, 2.5637077173107343, 2.5847217149936093, 2.6057357126764837, 2.6267497103593587, 2.6477637080422336, 2.6687777057251085, 2.6897917034079835, 2.7108057010908584, 2.731819698773733, 2.752833696456608, 2.7738476941394827, 2.7948616918223577, 2.8158756895052326, 2.8368896871881075, 2.857903684870982, 2.878917682553857, 2.8999316802367323, 2.9209456779196072, 2.941959675602482, 2.9629736732853567, 2.9839876709682316, 3.0050016686511065, 3.0260156663339814, 3.0470296640168564, 3.0680436616997313, 3.089057659382606, 3.1100716570654807, 3.1310856547483557, 3.1520996524312306, 3.1731136501141055, 3.1941276477969804, 3.215141645479855, 3.23615564316273, 3.257169640845605, 3.2781836385284797, 3.2991976362113546, 3.3202116338942296, 3.341225631577104, 3.362239629259979, 3.383253626942854, 3.404267624625729, 3.425281622308604, 3.4462956199914787, 3.467309617674353, 3.488323615357228, 3.509337613040103, 3.530351610722978, 3.551365608405853, 3.572379606088728, 3.5933936037716023, 3.6144076014544773, 3.635421599137352, 3.656435596820227, 3.6774495945031025, 3.6984635921859774, 3.7194775898688524, 3.7404915875517273, 3.761505585234602, 3.782519582917477, 3.8035335806003516, 3.8245475782832266, 3.8455615759661015, 3.8665755736489764, 3.8875895713318513, 3.9086035690147263, 3.9296175666976008, 3.9506315643804757, 3.9716455620633506, 3.9926595597462256, 4.0136735574291, 4.034687555111975, 4.05570155279485, 4.076715550477725, 4.0977295481606, 4.118743545843475, 4.13975754352635, 4.1607715412092245, 4.1817855388920995, 4.202799536574974, 4.223813534257849, 4.244827531940723, 4.265841529623598, 4.286855527306473, 4.307869524989348, 4.328883522672223, 4.349897520355098, 4.370911518037973, 4.391925515720848, 4.412939513403723, 4.433953511086598, 4.454967508769473, 4.475981506452348, 4.496995504135222, 4.518009501818097, 4.5390234995009715, 4.560037497183846, 4.581051494866721, 4.602065492549596, 4.623079490232471, 4.644093487915346, 4.665107485598221, 4.686121483281096, 4.707135480963971, 4.728149478646846, 4.749163476329721, 4.770177474012596, 4.791191471695471, 4.812205469378346, 5.085387439255719, 5.106401436938594, 5.127415434621469, 5.148429432304344, 5.1694434299872185, 5.190457427670093, 5.2114714253529675, 5.232485423035842, 5.253499420718717, 5.274513418401592, 5.295527416084467, 5.316541413767342, 5.337555411450217, 5.358569409133092, 5.379583406815967, 5.400597404498842, 5.421611402181717, 5.442625399864592]
y2 = [-0.17489399333249445, -0.27848476317808046, -0.29989085641034907, -0.2996597640653808, 0.29388777786073694, -0.29878959979806646, -0.29817364943935637, -0.2973993486197289, -0.29651570755810447, -0.2954948645869568, -0.2943444742201462, -0.29306068957991355, -0.2916518799826636, -0.29011066897845605, -0.2884439129085589, -0.2866493824696056, -0.28460542900102653, -0.2826571442192379, -0.2804903519011928, -0.27817906666651576, -0.27577674211587594, -0.27323811968110484, -0.27056512217364376, -0.2677892702459495, -0.26489747261560215, -0.2618718200808846, -0.2587376588373428, -0.2554981159687207, -0.2521391450273408, -0.24866174963401674, -0.24508080312868358, -0.2413976413637358, -0.23760237978229953, -0.23369646789354517, -0.2296920544945529, -0.22559017570734632, -0.22138848306858805, -0.2170823911558896, -0.21268194991746273, -0.20818876422701105, -0.20360805474672125, -0.19893582699659187, -0.19417108483784343, -0.1893207325844832, -0.18438932491680188, -0.17937736939572893, -0.17428833633150856, -0.16911976762813513, -0.16387369683243147, -0.15855530351041056, -0.15317120165644635, -0.1477151717170752, -0.14219842711624106, -0.1366187041292185, -0.13094505509333154, -0.12525417609398534, -0.11950895250633721, -0.11370737264111819, -0.10784907066406674, -0.10194331973143517, -0.09599047852480568, -0.08999336384324508, -0.08395731197791949, -0.07788561016187251, -0.07178093770980912, -0.06564565383210276, -0.05948202365672758, -0.053292365862272935, -0.04707912058790595, -0.0408448590834285, -0.034592261445046, -0.02832408339822046, -0.022043124778093854, -0.015752205230697335, -0.009454148027607474, -0.0031517699224263658, 0.0031521260655554715, 0.009454745941371794, 0.0157533039703916, 0.022045026254508184, 0.028327143732071965, 0.0345968850073433, 0.040851467850176605, 0.04708808950088521, 0.053303918157347856, 0.05949609004324046, 0.06566171826859984, 0.07179792005629962, 0.07790186596766321, 0.08397084695652597, 0.09000234418910072, 0.0959941078815303, 0.10194542234590104, 0.10785254834356962, 0.1137091560908198, 0.11951849476492515, 0.12527583497655048, 0.1309770306051448, 0.13661870412921767, 0.14219842711623945, 0.14771517171707563, 0.15317120165644713, 0.15855530351040992, 0.16387369683243236, 0.1691197676281353, 0.1742883363315086, 0.1793773693957285, 0.18438932491680163, 0.1893207325844828, 0.19417108483784407, 0.19893582699659274, 0.20360805474672186, 0.20818876422700902, 0.21268194991745992, 0.2170823911558873, 0.22142259007241988, 0.22561491871529946, 0.22970018994916164, 0.2337072004136354, 0.23762260189053788, 0.24141737075143566, 0.2450917289917178, 0.24867093997951825, 0.2521540262551526, 0.2555130592910183, 0.25874809478190286, 0.2618838410708835, 0.26490851610877497, 0.2678034243472967, 0.27057861409931167, 0.2732495815300367, 0.2757905486311947, 0.27820246053280523, 0.2805065648938426, 0.282679635494678, 0.28472129124907813, 0.2866893704412291, 0.28847514228785937, 0.29013320319423175, 0.2916576067941383, 0.29307731326830244, 0.2943506289142849, 0.2955059115462241, 0.2965239965251926, 0.29741051856427925, 0.2981773408868851, 0.2988093312042297, -0.2938877778607325, -0.28406137724102365, -0.27869012743522337, 0.27848476317808313, 0.28774542102218703, 0.1948843730146601, 0.2188485275818719, 0.24587776940822748, 0.2751972442698021, -0.2748687212283345, -0.24479604649179176, -0.2133881801153228, -0.1811858534388664, -0.14862698538550312, -0.11614962833688836, -0.08438506432738357, -0.05477275368794299, -0.032240262261240105, 0.05185381720490412, 0.07208388865104087, 0.10042240012629122, 0.13096602330320994, 0.16212995843793276, 0.19336385657292876, 0.22442256003846062, 0.2551742952552081, 0.2855373810481021]
z2 = [0.9008385399131996, -0.9380979647775209, -0.9233393381729095, -0.9210094205476991, 0.43604053877021626, 0.9134710877108979, 0.950139113772997, 0.49549071791407695, -0.9801405074904452, 0.9869080890723443, -0.9876724717138659, 0.987663072917269, -0.9452964646014579, 0.9620337397821722, 0.8844794075409858, 0.05013876371635637, -0.9811466334554276, -0.6824496027145411, -0.7163420460920841, -0.9868553366492926, -0.767345016796053, -0.8085824190174711, -0.9714090393186852, -0.7970836227880176, -0.9795456795222071, -0.8136358972463309, -0.8835172583558022, -0.8947853068305845, -0.9236716237136059, -0.8244736682212556, -0.7874317810157796, -0.9338669266150573, -0.9512680377393841, -0.7528157432848345, -0.7579389130442895, -0.8936393189098227, -0.9875912540064858, -0.880644066984375, -0.6401063182557362, -0.6988522470015288, -0.9405546324130515, -0.9875135009214225, -0.9252921690082477, -0.4434011868350803, -0.18449741494687535, -0.886829941933964, -0.9772485064317582, -0.9859909091383455, -0.9441266957664971, -0.4624501370966512, 0.9528750893887948, -0.7196121409962177, -0.9531491648632098, -0.9832071544338404, 0.9875661971662879, 0.9796562660377502, 0.9386654805123322, -0.013183205905960522, -0.9867518371972562, -0.8191971822195587, 0.8491621234399299, 0.9654820605458457, 0.9816144615990002, 0.9861184486467711, 0.9874908416631648, 0.9876805623656778, 0.9873316126468186, 0.986697599204872, 0.9858785676103959, 0.9848938166113397, 0.9836997245316758, 0.9821839734831609, 0.9801758100807315, 0.9777120819312816, 0.976868238482537, 0.984526812400976, 0.9845269307746456, 0.9768685754839224, 0.9777125805488167, 0.9801763672715285, 0.9821844312268214, 0.9836998641760942, 0.9848933645384209, 0.9858772229021093, 0.9866950944484474, 0.9873278140816761, 0.9876756013965708, 0.9874852440395594, 0.9861131938283625, 0.9816108777929695, 0.9654814699219721, 0.8491652345279307, -0.8191909614619228, -0.9867447647802569, -0.013178726694693046, 0.938664336348905, 0.9796490638569674, 0.9875561154667397, 0.9832071544338395, 0.9531491648632101, 0.719612140996233, -0.9528750893887912, 0.4624501370966147, 0.9441266957664991, 0.9859909091383456, 0.9772485064317596, 0.8868299419339634, 0.1844974149467401, 0.4434011868351206, 0.9252921690082508, 0.9875135009214235, 0.940554632413056, 0.6988522470014936, 0.6401063182557504, 0.8806440669843648, -0.9873910973636967, -0.893601478062046, -0.758130288070695, -0.7529400988252413, -0.9511609110527084, -0.9337690997262548, -0.7875554583013428, -0.8246041093047022, -0.923578789848213, -0.8947013829661918, -0.8836659831784237, -0.8137199674415136, -0.9793950736629681, -0.7971099867966105, -0.9716255210341596, -0.8084594327500895, -0.7672422720347151, -0.9872762771272058, -0.7161457125175263, -0.6820401687559688, -0.9869484021628372, 0.060114912689667906, 0.8842487512105972, 0.9619118430405836, -0.9452591585131124, 0.9876361931663517, -0.9876667159749235, 0.9869084238864297, -0.9801574101062112, 0.4963254584687391, 0.9500470276772067, 0.9134440585130736, -0.43604053877102456, -0.8697818945545679, -0.9327463875451452, 0.9380979647775277, 0.9455493391217722, -0.9354108566130968, -0.9555858491156347, -0.9681735644155438, -0.9764628775449539, 0.9912747897244912, 0.988841115056256, 0.9852511452682066, 0.979659624699124, 0.9702684089427748, 0.9527402039918781, 0.9146062435471375, 0.8105428353127574, 0.4473924787077477, 0.36412543420036075, 0.7906789861787451, 0.9105794656039878, 0.952489734272414, 0.9710648437848002, 0.9807408969261265, 0.9863710496432065, 0.9899173058797459, 0.9922872214419207];
















error()

















P = rand(10,7) 


p1 = Utils.PathConnect(P, 100, dim=1)[1]


@show size(p1)


p2 = Utils.PathConnect(transpose(P), 100,dim=2)[1]


@show size(p2)




@show isapprox(p1,transpose(p2))


#p1 .- transpose(p2) 

#@show P[3,:] p1[end,:] p2[:,end]




n = 123 
v = [1,0,3,1]

d =  Utils.PropDistributeBallsToBoxes(n, v)

@show d sum(d) n Utils.LA.normalize(d,1)*sum(v)


v = [1,1,3,1]

println.(Utils.PathConnect([-1.2,-0.6,0,0.6,1.2],10))

println() 

println.(Utils.PathConnect([-1.2,-0.6,0,0.6,1.2],10,v))




@show Utils.PathConnect(Float64[],100)
@show Utils.PathConnect(Float64[],100,v)




@show Utils.IdentifySectors([1,1,1,3,4,5,5,1]) 

@show Utils.IdentifyRanges([1,1,2,3,4,5,5,1])

@show Utils.IdentifyRanges([1,1,3,5,7,5,5,1,2,3])

@show Utils.IdentifyRanges(Utils.Random_Items(1:10))



for n in ["abc",1,10.5], d in [(),1,3,(1,2),(2,)]

#	println();	@show n d  

	Utils.nr2string(n,d)# |> println


end 



kPoints = Matrix{Float64}(undef, 0, 1)
NR_KPOINTS = 30

@show Utils.PathConnect(kPoints, NR_KPOINTS, dim=2)




using OrderedCollections: OrderedDict 

Utils.NT(Dict(:a=>2))#|>println

Utils.NT(OrderedDict(:a=>2))#|>println

Utils.NT((a=2,)) #|>println








#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


println()


include("utils_temp1.jl")

@show typeof(X1) 

obs = setdiff(keys(X1),["Energy"]) |> first 
@show obs 
println()

println()

@show size(X1[obs]["A"])
a=rand(3,5);b=rand(3,5);


@show cat(a,b,dims=2) |>size

@show Utils.catNewAx(2,a,b) |> size 
@show Utils.catNewAx(1,a,b) |> size 

@show Utils.catNewAx(2, a,Int.(round.(10b)), rand(3,5,1)) |> size 

println()
Energies = X1["Energy"]*range(1,-1,length=5)

@show Energies


function X(i)

	i==1 && return X1

	X2 = deepcopy(X1)

	X2["Energy"] = Energies[i]

#	for (k,v) in X2[obs]

	return X2 

end 


@show X(1)["Energy"]
@show X(5)["Energy"]


x = Utils.RecursiveMerge(X, 1:5; dim=2)[obs]
	
@show x["A"] |> size


@show X1[obs]["A"] |> size

#
#n = 1000
#
#vecs = map(1:100) do i 
#
#	L = rand(1:10) 
#
#	s = ones(Int,L)
#
#	s[rand(1:L)] = n
#
#	return rand(s...)
#
#end 
#
#Utils.VecsToMat(vecs...; dim=2);
#@time Utils.VecsToMat(vecs...; dim=2);
#@btime Utils.VecsToMat(vecs...; dim=2); 
#



d1 = Dict(1=>2,:a=>:b,3=>rand(2),"s"=>rand(1,1))

d2 = [Dict(1=>3)]

println(Utils.dict_keepkeys(d1,d2...))


@show Utils.uniqlogsp(2,6,6,3; Trunc=true)
#@show Utils.uniqlogsp2(2,6,6,3; Trunc=true)
@show Utils.uniqlinsp(2,6,6,3; Trunc=true)

result = [1, 2, 4, 6, 7] 

for x in [(5,), (rand(5,3),1), (1:5,), ([1,2,3,4,5],)]

	for y in [(7,), (rand(5,7),2), (1:7,), ([1,2,7],)]

#		println("First: ", typeof.(x)...)
#		println("Second: ", typeof.(y)...)
		
		@assert result==Utils.RescaleInds(x..., y...) 
	
#		println() 
	end 
end 



#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


println()
println("----------------")
println("DistributeBallsToBoxes")
println()






for item in Utils.DistributeBallsToBoxes(1,3)

	println(item,"\n")

end 

g14(x) = [x,-x]
f = g14
#Utils.inv_abs 


@time Utils.DistributeBallsToBoxes(2,3,f)
@time Utils.DistributeBallsToBoxes(2,3,f)

@time Utils.DistributeBallsToBoxes(2,3,g14)
@time Utils.DistributeBallsToBoxes(2,3,g14) 

for item in Utils.DistributeBallsToBoxes(2,3,f)

#	println(item," ",sum(abs,item),"\n")
	@assert sum(abs,item) == 2

end 

println("----------------")
println("\nCombsOfVecs\n")

@testset "CombsOfVecs" begin 


vecs= rand(2,5)
coeffs = rand(5,3)


cv = Utils.CombsOfVecs(vecs, coeffs; dim=2)


@test size(cv)==(2,3)



cv2 = Utils.CombsOfVecs(vecs[:,2], coeffs[1:1,:]; dim=2)

@test size(cv2)==(2,3)



cv3 = Utils.CombsOfVecs(vecs, coeffs[:,1:1]; dim=2)


@test size(cv3)==(2,1)




cv4 = Utils.CombsOfVecs(vecs, coeffs[:,1]; dim=2)


@test size(cv4)==(2,) && isapprox(cv4,cv3[:])






end







