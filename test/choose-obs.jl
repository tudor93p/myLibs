import myLibs: ReadWrite, ComputeTasks 
import JLD 

fn = x->"/media/tudor/Tudor/Work/2020_Snake-states/SnakeStates/Data/DeviceWithLeads/GreensFcts/075-37-2p000-0p030-2-0p005-0p00-0p400-0p400-0p002/AB/18-1-1p0-0p0-m1-1p00-m1p0-1p0/18-1-1p0-0p0-1-1p00-1p0-m1p0/$x"

obs = "QP-SiteVectorTransmission"


@show isfile(fn(obs)*".jld")

@time Data = ReadWrite.Read_NamesVals(fn, obs, "jld")



@show typeof(Data) typeof(Data[obs]) 

for P in [Dict("obs_i"=>1),:obs_i=>1.3,1,Dict(),7,1.0,"1.12342",(obs_i=1.1,x=2)]

	@assert 1==ComputeTasks.parse_obs_i(P, rand(5), "first")

	@time ComputeTasks.choose_obs_i(Data[obs]; P=P, f="first")
	@time ComputeTasks.choose_obs_i(Data[obs]; P=P, f="first")


end 
